
#include <stdio.h>
#include "cow.h"
#if (COW_MPI)
#include <mpi.h>
#endif
#define KILOBYTES (1<<10)
#define MEGABYTES (1<<20)

static void divergence(double *result, double **args, int **s, cow_domain *d)
{
#define M(i,j,k) ((i)*s[0][0] + (j)*s[0][1] + (k)*s[0][2])
  double *fx = &args[0][0];
  double *fy = &args[0][1];
  double *fz = &args[0][2];
  *result = ((fx[M(1,0,0)] + fx[M(1,1,0)] + fx[M(1,0,1)] + fx[M(1,1,1)]) -
	     (fx[M(0,0,0)] + fx[M(0,1,0)] + fx[M(0,0,1)] + fx[M(0,1,1)])) / 4.0
    +       ((fy[M(0,1,0)] + fy[M(0,1,1)] + fy[M(1,1,0)] + fy[M(1,1,1)]) -
	     (fy[M(0,0,0)] + fy[M(0,0,1)] + fy[M(1,0,0)] + fy[M(1,0,1)])) / 4.0
    +       ((fz[M(0,0,1)] + fz[M(1,0,1)] + fz[M(0,1,1)] + fz[M(1,1,1)]) -
	     (fz[M(0,0,0)] + fz[M(1,0,0)] + fz[M(0,1,0)] + fz[M(1,1,0)])) / 4.0;
#undef M
}
static void curl(double *result, double **args, int **s, cow_domain *d)
{
// http://en.wikipedia.org/wiki/Five-point_stencil
#define diff5(f,s) ((-f[2*s] + 8*f[s] - 8*f[-s] + f[-2*s]) / 12.0)
  double *f0 = &args[0][0];
  double *f1 = &args[0][1];
  double *f2 = &args[0][2];
  result[0] = diff5(f2, s[0][1]) - diff5(f1, s[0][2]);
  result[1] = diff5(f0, s[0][2]) - diff5(f2, s[0][0]);
  result[2] = diff5(f1, s[0][0]) - diff5(f0, s[0][1]);
#undef diff5
}

int main(int argc, char **argv)
{
#if (COW_MPI)
  {
    int rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank != 0) freopen("/dev/null", "w", stdout);
    printf("was compiled with MPI support\n");
  }
#endif

  if (argc == 3) {
    printf("running on input file %s\n", argv[1]);
  }
  else {
    printf("usage: $> milos infile.h5 outfile.h5\n");
    goto done;
  }

  cow_domain *domain = cow_domain_new();
  cow_dfield *vel = cow_dfield_new(domain, "prim");
  cow_dfield *mag = cow_dfield_new(domain, "prim");
  cow_dfield *divB = cow_dfield_new(domain, "divB");
  cow_dfield *divV = cow_dfield_new(domain, "divV");
  cow_dfield *curlB = cow_dfield_new(domain, "curlB");
  cow_dfield *curlV = cow_dfield_new(domain, "curlV");

  int collective = getenv("COW_HDF5_COLLECTIVE") ? atoi(getenv("COW_HDF5_COLLECTIVE")) : 0;
  printf("COW_HDF5_COLLECTIVE: %d\n", collective);

  cow_domain_readsize(domain, argv[1], "prim/vx");
  cow_domain_setguard(domain, 2);
  cow_domain_commit(domain);

  cow_domain_setchunk(domain, 1);
  cow_domain_setcollective(domain, collective);
  cow_domain_setalign(domain, 4*KILOBYTES, 4*MEGABYTES);

  cow_dfield_addmember(vel, "vx");
  cow_dfield_addmember(vel, "vy");
  cow_dfield_addmember(vel, "vz");
  cow_dfield_addmember(mag, "Bx");
  cow_dfield_addmember(mag, "By");
  cow_dfield_addmember(mag, "Bz");

  cow_dfield_addmember(divB, "divB");
  cow_dfield_addmember(divV, "divV");
  cow_dfield_addmember(curlB, "fx");
  cow_dfield_addmember(curlB, "fy");
  cow_dfield_addmember(curlB, "fz");
  cow_dfield_addmember(curlV, "fx");
  cow_dfield_addmember(curlV, "fy");
  cow_dfield_addmember(curlV, "fz");

  cow_dfield_commit(vel);
  cow_dfield_commit(mag);
  cow_dfield_commit(divB);
  cow_dfield_commit(divV);
  cow_dfield_commit(curlB);
  cow_dfield_commit(curlV);

  cow_dfield_read(vel, argv[1]);
  cow_dfield_read(mag, argv[1]);

  cow_dfield_transform(divB, &mag, 1, divergence);
  cow_dfield_transform(divV, &vel, 1, divergence);
  cow_dfield_transform(curlB, &mag, 1, curl);
  cow_dfield_transform(curlV, &vel, 1, curl);

  cow_dfield_write(divB, argv[2]);
  cow_dfield_write(divV, argv[2]);
  cow_dfield_write(curlB, argv[2]);
  cow_dfield_write(curlV, argv[2]);

  cow_dfield_del(vel);
  cow_dfield_del(mag);
  cow_dfield_del(divB);
  cow_dfield_del(divV);
  cow_dfield_del(curlB);
  cow_dfield_del(curlV);
  cow_domain_del(domain);

 done:
#if (COW_MPI)
  MPI_Finalize();
#endif
  return 0;
}
