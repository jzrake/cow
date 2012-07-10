
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#define COW_PRIVATE_DEFS
#include "cow.h"

#define MODULE "sampling"

int cow_domain_indexatposition(cow_domain *d, int dim, double x)
// -----------------------------------------------------------------------------
// d: 0,1,2 for x,y,z
//
// r is a global position, i.e. not relative to this subgrid
// The return value is the integer index, which is relative to this subgrid
//
// Within the zone i+Ng, the value (x-x0)/dx ranges from i to i+1.
// Then we correct for ghosts by adding Ng.
// -----------------------------------------------------------------------------
{
  if (dim >= 3 || !d->committed) return 0.0;
  return d->n_ghst + (int)((x - d->loc_lower[dim]) / d->dx[dim]);
}

double cow_domain_positionatindex(cow_domain *d, int dim, int index)
// -----------------------------------------------------------------------------
// d: 0,1,2 for x,y,z
//
// r is a global position, i.e. not relative to this subgrid
// The return value is the integer index, which is relative to this subgrid
//
// Within the zone i+Ng, the value (x-x0)/dx ranges from i to i+1.
// Then we correct for ghosts by adding Ng.
// -----------------------------------------------------------------------------
{
  if (dim >= 3 || !d->committed) return 0.0;
  return d->loc_lower[dim] + d->dx[dim] * (index - d->n_ghst + 0.5);
}

void cow_dfield_sample(cow_dfield *f, double *x, double *P)
{
#define M(i,j,k) ((i)*s[0] + (j)*s[1] + (k)*s[2])
  cow_domain *d = f->domain;
  int *s = f->stride;
  int i = cow_domain_indexatposition(d, 0, x[0]);
  int j = cow_domain_indexatposition(d, 1, x[1]);
  int k = cow_domain_indexatposition(d, 2, x[2]);
  double x0 = cow_domain_positionatindex(d, 0, i-1);
  double y0 = cow_domain_positionatindex(d, 1, j-1);
  double z0 = cow_domain_positionatindex(d, 2, k-1);
  double *A = (double*) f->data;
  double *P000 = &A[M(i-1,j-1,k-1)];
  double *P001 = &A[M(i-1,j-1,k+1)];
  double *P010 = &A[M(i-1,j+1,k-1)];
  double *P011 = &A[M(i-1,j+1,k+1)];
  double *P100 = &A[M(i+1,j-1,k-1)];
  double *P101 = &A[M(i+1,j-1,k+1)];
  double *P110 = &A[M(i+1,j+1,k-1)];
  double *P111 = &A[M(i+1,j+1,k+1)];
  double delx[3] = {
    0.5 * (x[0] - x0) / d->dx[0],
    0.5 * (x[1] - y0) / d->dx[1],
    0.5 * (x[2] - z0) / d->dx[2] };
  // ---------------------------------------------------------------------------
  // See http://en.wikipedia.org/wiki/Trilinear_interpolation
  // ---------------------------------------------------------------------------
  for (int q=0; q<f->n_members; ++q) {
    double i1 = P000[q] * (1.0 - delx[2]) + P001[q] * delx[2];
    double i2 = P010[q] * (1.0 - delx[2]) + P011[q] * delx[2];
    double j1 = P100[q] * (1.0 - delx[2]) + P101[q] * delx[2];
    double j2 = P110[q] * (1.0 - delx[2]) + P111[q] * delx[2];
    double w1 = i1 * (1.0 - delx[1]) + i2 * delx[1];
    double w2 = j1 * (1.0 - delx[1]) + j2 * delx[1];
    P[q] = w1 * (1.0 - delx[0]) + w2 * delx[0];
  }
#undef M
}
