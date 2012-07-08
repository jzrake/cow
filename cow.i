/* -*- C -*-  (not really, but good for syntax highlighting) */

%module cow
%{
#define SWIG_FILE_WITH_INIT
#include "cow.h"
  void test_trans(double *result, double **args, int **strides, void *udata);
  void getarray3(cow_dfield *f, double **data, int *dim1, int *dim2, int *dim3)
  {
    *dim1 = cow_domain_getnumlocalzonesincguard(cow_dfield_getdomain(f), 0);
    *dim2 = cow_domain_getnumlocalzonesincguard(cow_dfield_getdomain(f), 1);
    *dim3 = cow_domain_getnumlocalzonesincguard(cow_dfield_getdomain(f), 2);
    *data = (double*) cow_dfield_getdata(f);
  }
  %}

%typemap(in) (const char *name)
{
  $1 = PyString_AsString($input);
}

%include "numpy.i"
%typemap(in,numinputs=0)
(cow_dfield *f, DATA_TYPE** ARGOUTVIEW_ARRAY3,
 DIM_TYPE* DIM1, DIM_TYPE* DIM2, DIM_TYPE* DIM3)
(cow_dfield *f, DATA_TYPE* data_temp,
 DIM_TYPE dim1_temp, DIM_TYPE dim2_temp, DIM_TYPE dim3_temp)
{
  $1 = &f;
  $2 = &data_temp;
  $3 = &dim1_temp;
  $4 = &dim2_temp;
  $5 = &dim3_temp;
}
%typemap(argout,
         fragment="NumPy_Backward_Compatibility")
  (cow_dfield *f, DATA_TYPE** ARGOUTVIEW_ARRAY3,
   DIM_TYPE* DIM1, DIM_TYPE* DIM2, DIM_TYPE* DIM3)
{
  npy_intp dims[3] = { *$3, *$4, *$5 };
  PyObject * array = PyArray_SimpleNewFromData(3, dims, DATA_TYPECODE, (void*)(*$2));
  if (!array) SWIG_fail;
  $result = SWIG_Python_AppendOutput($result,array);
}

%init %{
import_array();
%}

%apply (double** ARGOUTVIEW_ARRAY3, int* DIM1, int* DIM2, int* DIM3)
{
  (double **data, int *dim1, int *dim2, int *dim3) };
%apply (cow_dfield *f, double** ARGOUTVIEW_ARRAY3, int* DIM1, int* DIM2, int* DIM3)
{
  (cow_dfield *f, double **data, int *dim1, int *dim2, int *dim3) };

%include "cow.h"
void getarray3(cow_dfield *f, double **data, int *dim1, int *dim2, int *dim3);

%callback("%(upper)s");
void test_trans(double *result, double **args, int **strides, void *udata);
%nocallback;
