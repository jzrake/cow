/* -*- C -*-  (not really, but good for syntax highlighting) */

%module cow
%{
#define SWIG_FILE_WITH_INIT
#include <numpy/arrayobject.h>
#include "cow.h"
  void setarray1(cow_dfield *f, void *x, int n0, int n1)
  {
    cow_dfield_setbuffer(f, x);
  }
  void setarray2(cow_dfield *f, void *x, int n0, int n1, int n2)
  {
    cow_dfield_setbuffer(f, x);
  }
  void setarray3(cow_dfield *f, void *x, int n0, int n1, int n2, int n3)
  {
    cow_dfield_setbuffer(f, x);
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
%apply(double *IN_ARRAY3, int DIM1, int DIM2)
{(void *x, int n0, int n1)};
%apply(double *IN_ARRAY3, int DIM1, int DIM2, int DIM3)
{(void *x, int n0, int n1, int n2)};
%apply(double *IN_ARRAY4, int DIM1, int DIM2, int DIM3, int DIM4)
{(void *x, int n0, int n1, int n2, int n3)};

%include "cow.h"
extern void setarray1(cow_dfield *f, void *x, int n0, int n1);
extern void setarray2(cow_dfield *f, void *x, int n0, int n1, int n2);
extern void setarray3(cow_dfield *f, void *x, int n0, int n1, int n2, int n3);

%clear(double x[3]);
%clear(double *x, int n0, int n1);
%clear(double *x, int n0, int n1, int n2);
%clear(double *x, int n0, int n1, int n2, int n3);

