/* -*- C -*-  (not really, but good for syntax highlighting) */

// http://www.scipy.org/Cookbook/SWIG_NumPy_examples

%module cow
%{
#define SWIG_FILE_WITH_INIT

  // Swig source is littered with 'no prototype' warnings on icc:
#ifdef __INTEL_COMPILER
#pragma warning disable 1418
#endif

#include <numpy/arrayobject.h>
#include "cow.h"
void setarray1(cow_dfield *f, double *x, int n0, int n1)
{
  cow_dfield_setbuffer(f, x);
}
void setarray2(cow_dfield *f, double *x, int n0, int n1, int n2)
{
  cow_dfield_setbuffer(f, x);
}
void setarray3(cow_dfield *f, double *x, int n0, int n1, int n2, int n3)
{
  cow_dfield_setbuffer(f, x);
}

void cow_trans_divcorner(double *result, double **args, int **s, void *u);
void cow_trans_div5(double *result, double **args, int **s, void *u);
void cow_trans_rot5(double *result, double **args, int **s, void *u);
void cow_trans_component(double *result, double **args, int **s, void *u);

  %}

%typemap(in) (const char *name)
{
  $1 = PyString_AsString($input);
}
%init %{
  import_array();
  %}


%include "numpy.i"

 //%apply(double IN_ARRAY1[ANY]){(double x[3])};
%apply(double ARGOUT_ARRAY1[ANY]){(double x[3])};

%apply(double *IN_ARRAY2, int DIM1, int DIM2)
{(double *x, int n0, int n1)};
%apply(double *IN_ARRAY3, int DIM1, int DIM2, int DIM3)
{(double *x, int n0, int n1, int n2)};
%apply(double *IN_ARRAY4, int DIM1, int DIM2, int DIM3, int DIM4)
{(double *x, int n0, int n1, int n2, int n3)};
%apply(double **ARGOUTVIEW_ARRAY1, int *DIM1)
{(double **x, int *n0)};
%apply(double **ARGOUTVIEW_ARRAY2, int *DIM1, int *DIM2)
{(double **x, int *n0, int *n1)};

%include "cow.h"

extern void setarray1(cow_dfield *f, double *x, int n0, int n1);
extern void setarray2(cow_dfield *f, double *x, int n0, int n1, int n2);
extern void setarray3(cow_dfield *f, double *x, int n0, int n1, int n2, int n3);

%clear(double *x, int n0, int n1);
%clear(double *x, int n0, int n1, int n2);
%clear(double *x, int n0, int n1, int n2, int n3);
%clear(double **x, int *n0);
%clear(double **x, int *n1, int *n0);

%constant void cow_trans_divcorner(double *result, double **args, int **s, void *u);
%constant void cow_trans_div5(double *result, double **args, int **s, void *u);
%constant void cow_trans_rot5(double *result, double **args, int **s, void *u);
%constant void cow_trans_component(double *result, double **args, int **s, void *u);
