/* -*- C -*-  (not really, but good for syntax highlighting) */

%module cow
%{
#define SWIG_FILE_WITH_INIT
#include <numpy/arrayobject.h>
#include "cow.h"
  void testfunc1(cow_domain *d, double x[3])
  {
    printf("testfunc1: %f %f %f\n", x[0], x[1], x[2]);
  }
  void testfunc2(cow_domain *d, double *x, int n0, int n1, int n2)
  {
    printf("testfunc2: %d %d %d\n", n0, n1, n2);
  }
  %}

%typemap(in) (const char *name)
{
  $1 = PyString_AsString($input);
}
%init %{
  import_array();
  %}


%include "numpy.i"

%apply(double IN_ARRAY1[ANY]){(double x[3])};
%apply(double *IN_ARRAY3, int DIM1, int DIM2, int DIM3){(double *x, int n0, int n1, int n2)};
%include "cow.h"
%clear(double x[3]);
%clear(double *x, int n0, int n1, int n2);

