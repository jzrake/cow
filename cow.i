/* -*- C -*-  (not really, but good for syntax highlighting) */

%module cow
%{
#define SWIG_FILE_WITH_INIT
#include <numpy/arrayobject.h>
#include "cow.h"
  void test_trans(double *result, double **args, int **strides, void *udata);
  %}

%typemap(in) (const char *name)
{
  $1 = PyString_AsString($input);
}
%init %{
  import_array();
  %}

%include "cow.h"
%native(cow_dfield_getarray) PyObject *_wrap_cow_dfield_getarray(PyObject *SWIGUNUSEDPARM(self), PyObject *args);

%callback("%(upper)s");
void test_trans(double *result, double **args, int **strides, void *udata);
%nocallback;








%{
SWIGINTERN PyObject *_wrap_cow_dfield_getarray(PyObject *SWIGUNUSEDPARM(self), PyObject *args)
{
  PyObject *resultobj = 0;
  void *argp1 = 0;
  int res1 = 0;
  PyObject *obj0 = 0;

  if (!PyArg_ParseTuple(args, (char*)"O:cow_dfield_getarray",&obj0)) SWIG_fail;
  res1 = SWIG_ConvertPtr(obj0, &argp1, SWIGTYPE_p_cow_dfield, 0 | 0);

  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1),
                        "in method '" "cow_dfield_getarray" "', argument " "1"
                        " of type '" "cow_dfield *""'");
  }

  cow_dfield *f = (cow_dfield *)(argp1);
  cow_domain *d = cow_dfield_getdomain(f);
  double *data = (double*) cow_dfield_getbuffer(f);
  npy_intp dims[4];
  int ndims = cow_domain_getndim(d);
  dims[0] = cow_domain_getnumlocalzonesincguard(d, 0);
  dims[1] = cow_domain_getnumlocalzonesincguard(d, 1);
  dims[2] = cow_domain_getnumlocalzonesincguard(d, 2);
  dims[ndims] = cow_dfield_getnmembers(f);

  resultobj = SWIG_Py_Void();
  PyObject *array = PyArray_SimpleNewFromData(ndims+1, dims, NPY_DOUBLE, (void*)data);
  if (!array) SWIG_fail;
  resultobj = SWIG_Python_AppendOutput(resultobj, array);
  return resultobj;
 fail:
  return NULL;
}
 %}
