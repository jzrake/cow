
/*------------------------------------------------------------------------------
 * FILE: lua_fft.cpp
 *
 * AUTHOR: Jonathan Zrake, NYU CCPP: zrake@nyu.edu
 *
 * REFERENCES:
 *
 * http://www.cs.sandia.gov/~sjplimp/docs/fft/README.html
 *
 * NOTES:
 *
 * This code wraps the parallel FFT and remapping routines of Steve Plimpton at
 * Sandia National Labs.
 *
 * - The option to 'SCALED_YES' to divide by N is only regarded for forward
 *   transforms.
 *
 * - The order of indices provided to the FFT's is fast,mid,slow varying. For C
 *   arrays, this means it gets called with Nz, Ny, Nx.
 *
 *------------------------------------------------------------------------------
 */


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#define COW_PRIVATE_DEFS
#include "cow.h"

#include "fft_3d.h"
#define MODULE "fft"
#define NBINS 128




// Private functions used in this module
// -----------------------------------------------------------------------------
#define SCALED_NOT 0
#define SCALED_YES 1

#define PERMUTE_NONE 0
#define FFT_FWD (+1)
#define FFT_REV (-1)

static struct fft_plan_3d *call_fft_plan_3d(cow_domain *d, int *nbuf);
static double k_at(cow_domain *d, int i, int j, int k, double *khat);
static double khat_at(cow_domain *d, int i, int j, int k, double *khat);
static double cnorm(FFT_DATA z);

void cow_fft_pspecvecfield(cow_dfield *f, const char *fout, const char *gout)
{
  if (!f->committed) return;
  if (f->n_members != 3) {
    printf("[%s] error: need a 3-component field for pspecvectorfield", MODULE);
    return;
  }

  clock_t start = clock();
  int nx = cow_domain_getnumlocalzonesinterior(f->domain, 0);
  int ny = cow_domain_getnumlocalzonesinterior(f->domain, 1);
  int nz = cow_domain_getnumlocalzonesinterior(f->domain, 2);
  int Nx = cow_domain_getnumglobalzones(f->domain, 0);
  int Ny = cow_domain_getnumglobalzones(f->domain, 1);
  int Nz = cow_domain_getnumglobalzones(f->domain, 2);
  int ng = cow_domain_getguard(f->domain);
  int nbuf;
  int ntot = nx * ny * nz;
  int I0[3] = { ng, ng, ng };
  int I1[3] = { nx + ng, ny + ng, nz + ng };

  double *input = (double*) malloc(3 * ntot * sizeof(double));
  cow_dfield_extract(f, I0, I1, input);

  struct fft_plan_3d *plan = call_fft_plan_3d(f->domain, &nbuf);
  FFT_DATA *fx = (FFT_DATA*) malloc(nbuf * sizeof(FFT_DATA));
  FFT_DATA *fy = (FFT_DATA*) malloc(nbuf * sizeof(FFT_DATA));
  FFT_DATA *fz = (FFT_DATA*) malloc(nbuf * sizeof(FFT_DATA));
  printf("[%s] need to allocate %d zones\n", MODULE, nbuf);
  for (int i=0; i<nbuf; ++i) {
    fx[i].re = input[3*i + 0];
    fy[i].re = input[3*i + 1];
    fz[i].re = input[3*i + 2];
    fx[i].im = 0.0;
    fy[i].im = 0.0;
    fz[i].im = 0.0;
  }
  free(input);
  FFT_DATA *gx = (FFT_DATA*) malloc(nbuf * sizeof(FFT_DATA));
  FFT_DATA *gy = (FFT_DATA*) malloc(nbuf * sizeof(FFT_DATA));
  FFT_DATA *gz = (FFT_DATA*) malloc(nbuf * sizeof(FFT_DATA));
  fft_3d(fx, gx, FFT_FWD, plan);
  fft_3d(fy, gy, FFT_FWD, plan);
  fft_3d(fz, gz, FFT_FWD, plan);
  free(fx);
  free(fy);
  free(fz);
  cow_histogram *hist = cow_histogram_new();
  cow_histogram_setlower(hist, 0, 1.0);
  cow_histogram_setupper(hist, 0, 0.5*sqrt(Nx*Nx + Ny*Ny + Nz*Nz));
  cow_histogram_setnbins(hist, 0, NBINS);
  cow_histogram_setbinmode(hist, COW_HIST_BINMODE_DENSITY);
  cow_histogram_setspacing(hist, COW_HIST_SPACING_LOG);
  cow_histogram_setdomaincomm(hist, f->domain);
  cow_histogram_commit(hist);
  cow_histogram_setnickname(hist, cow_dfield_getname(f));
  for (int i=0; i<nx; ++i) {
    for (int j=0; j<ny; ++j) {
      for (int k=0; k<nz; ++k) {
        const int m = i*ny*nz + j*nz + k;
	double kvec[3];
        double khat[3];
	khat_at(f->domain, i, j, k, khat);
	// ---------------------------------------------------------------------
	// Here we are taking the complex norm (absolute value squared) of the
	// vector-valued Fourier amplitude corresponding to the wave-vector, k.
	//
	//                        P(k) = |\vec{f}_\vec{k}|^2
	//
	// ---------------------------------------------------------------------
	const double Kijk = k_at(f->domain, i, j, k, kvec);
	const double Pijk = cnorm(gx[m]) + cnorm(gy[m]) + cnorm(gz[m]);
	cow_histogram_addsample1(hist, Kijk, Pijk);
      }
    }
  }
  free(gx);
  free(gy);
  free(gz);
  cow_histogram_dumphdf5(hist, fout, gout);
  cow_histogram_del(hist);
  fft_3d_destroy_plan(plan);
  printf("[%s] %s took %3.2f seconds\n",
	 MODULE, __FUNCTION__, (double) (clock() - start) / CLOCKS_PER_SEC);
}


struct fft_plan_3d *call_fft_plan_3d(cow_domain *d, int *nbuf)
{
#if (COW_MPI)
  const int i0 = cow_domain_getglobalstartindex(d, 0);
  const int i1 = cow_domain_getnumlocalzonesinterior(d, 0) + i0 - 1;
  const int j0 = cow_domain_getglobalstartindex(d, 1);
  const int j1 = cow_domain_getnumlocalzonesinterior(d, 1) + j0 - 1;
  const int k0 = cow_domain_getglobalstartindex(d, 2);
  const int k1 = cow_domain_getnumlocalzonesinterior(d, 2) + k0 - 1;
  const int Nx = cow_domain_getnumglobalzones(d, 0);
  const int Ny = cow_domain_getnumglobalzones(d, 1);
  const int Nz = cow_domain_getnumglobalzones(d, 2);
  return fft_3d_create_plan(d->mpi_cart,
                            Nz, Ny, Nx,
                            k0,k1, j0,j1, i0,i1,
                            k0,k1, j0,j1, i0,i1,
                            SCALED_YES, PERMUTE_NONE, nbuf);
#else
  return NULL;
#endif
}

double k_at(cow_domain *d, int i, int j, int k, double *kvec)
// -----------------------------------------------------------------------------
// Here, we populate the wave vectors on the Fourier lattice. The convention
// used by FFTW is the same as that used by numpy, described at the link
// below. For N odd, the (positive) Nyquist frequency is placed in the middle
// bin.
//
// http://docs.scipy.org/doc/numpy/reference/generated/numpy.fft.fftfreq.html
// -----------------------------------------------------------------------------
{
  const int Nx = cow_domain_getnumglobalzones(d, 0);
  const int Ny = cow_domain_getnumglobalzones(d, 1);
  const int Nz = cow_domain_getnumglobalzones(d, 2);
  i += cow_domain_getglobalstartindex(d, 0);
  j += cow_domain_getglobalstartindex(d, 1);
  k += cow_domain_getglobalstartindex(d, 2);
  kvec[0] = (Nx % 2 == 0) ?
    ((i<  Nx   /2) ? i : i-Nx):  // N even
    ((i<=(Nx-1)/2) ? i : i-Nx);  // N odd
  kvec[1] = (Ny % 2 == 0) ?
    ((j<  Ny   /2) ? j : j-Ny):
    ((j<=(Ny-1)/2) ? j : j-Ny);
  kvec[2] = (Nz % 2 == 0) ?
    ((k<  Nz   /2) ? k : k-Nz):
    ((k<=(Nz-1)/2) ? k : k-Nz);
  return sqrt(kvec[0]*kvec[0] + kvec[1]*kvec[1] + kvec[2]*kvec[2]);
}
double khat_at(cow_domain *d, int i, int j, int k, double *khat)
{
  const double k0 = k_at(d, i, j, k, khat);
  if (fabs(k0) > 1e-12) { // don't divide by zero
    khat[0] /= k0;
    khat[1] /= k0;
    khat[2] /= k0;
  }
  return k0;
}
double cnorm(FFT_DATA z)
// http://www.cplusplus.com/reference/std/complex/norm
{
  return z.re*z.re + z.im*z.im;
}
