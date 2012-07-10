
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#define COW_PRIVATE_DEFS
#include "cow.h"

#define MODULE "sampling"

static void _sample1(cow_dfield *f, double *x, double *P);
static void _sample2(cow_dfield *f, double *x, double *P);
static void _sample3(cow_dfield *f, double *x, double *P);


void cow_dfield_sample(cow_dfield *f, double *x, double *P)
{
  switch (f->domain->n_dims) {
  case 1: _sample1(f, x, P); break;
  case 2: _sample2(f, x, P); break;
  case 3: _sample3(f, x, P); break;
  default: break;
  }
}
void _sample1(cow_dfield *f, double *x, double *P)
{
#define M(i) ((i)*s[0])
  cow_domain *d = f->domain;
  int *s = f->stride;
  int i = cow_domain_indexatposition(d, 0, x[0]);
  double x0 = cow_domain_positionatindex(d, 0, i-1);
  double *A = (double*) f->data;
  double *P0 = &A[M(i-1)];
  double *P1 = &A[M(i+1)];
  double delx[1] = { 0.5 * (x[0] - x0) / d->dx[0] };
  for (int q=0; q<f->n_members; ++q) {
    double b1 = P0[q];
    double b2 = P1[q] - P0[q];
    P[q] = b1 + b2*delx[0];
  }
#undef M
}
void _sample2(cow_dfield *f, double *x, double *P)
{
#define M(i,j) ((i)*s[0] + (j)*s[1])
  cow_domain *d = f->domain;
  int *s = f->stride;
  int i = cow_domain_indexatposition(d, 0, x[0]);
  int j = cow_domain_indexatposition(d, 1, x[1]);
  double x0 = cow_domain_positionatindex(d, 0, i-1);
  double y0 = cow_domain_positionatindex(d, 1, j-1);
  double *A = (double*) f->data;
  double *P00 = &A[M(i-1,j-1)];
  double *P01 = &A[M(i-1,j+1)];
  double *P10 = &A[M(i+1,j-1)];
  double *P11 = &A[M(i+1,j+1)];
  double delx[2] = {
    0.5 * (x[0] - x0) / d->dx[0],
    0.5 * (x[1] - y0) / d->dx[1] };
  for (int q=0; q<f->n_members; ++q) {
    double b1 = P00[q];
    double b2 = P10[q] - P00[q];
    double b3 = P01[q] - P00[q];
    double b4 = P00[q] - P10[q] - P01[q] + P11[q];
    P[q] = b1 + b2*delx[0] + b3*delx[1] + b4*delx[0]*delx[1];
  }
#undef M
}
void _sample3(cow_dfield *f, double *x, double *P)
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
