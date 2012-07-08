/* -*- C -*-  (not really, but good for syntax highlighting) */

%module cow
%{
#define SWIG_FILE_WITH_INIT
#include <numpy/arrayobject.h>
#include "cow.h"
  void test_trans(double *result, double **args, int **strides, void *udata);
  void getarray3(cow_dfield *f, double **data, int *dim1, int *dim2, int *dim3)
  {
    *dim1 = cow_domain_getnumlocalzonesincguard(cow_dfield_getdomain(f), 0);
    *dim2 = cow_domain_getnumlocalzonesincguard(cow_dfield_getdomain(f), 1);
    *dim3 = cow_domain_getnumlocalzonesincguard(cow_dfield_getdomain(f), 2);
    *data = (double*) cow_dfield_getdata(f);
  }
  SWIGINTERN PyObject *_wrap_getarray3(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
    PyObject *resultobj = 0;
    cow_dfield *arg1 = (cow_dfield *) 0 ;
    double **arg2 = (double **) 0 ;
    int *arg3 = (int *) 0 ;
    int *arg4 = (int *) 0 ;
    int *arg5 = (int *) 0 ;
    void *argp1 = 0 ;
    int res1 = 0 ;
    double *data_temp2 ;
    int dim1_temp2 ;
    int dim2_temp2 ;
    int dim3_temp2 ;
    PyObject * obj0 = 0 ;

    {
      arg2 = &data_temp2;
      arg3 = &dim1_temp2;
      arg4 = &dim2_temp2;
      arg5 = &dim3_temp2;
    }
    if (!PyArg_ParseTuple(args,(char *)"O:getarray3",&obj0)) SWIG_fail;
    res1 = SWIG_ConvertPtr(obj0, &argp1,SWIGTYPE_p_cow_dfield, 0 |  0 );
    if (!SWIG_IsOK(res1)) {
      SWIG_exception_fail(SWIG_ArgError(res1), "in method '" "getarray3" "', argument " "1"" of type '" "cow_dfield *""'");
    }
    arg1 = (cow_dfield *)(argp1);
    getarray3(arg1,arg2,arg3,arg4,arg5);
    resultobj = SWIG_Py_Void();
    {
      npy_intp dims[3] = {
        *arg3, *arg4, *arg5
      };
      PyObject * array = PyArray_SimpleNewFromData(3, dims, NPY_DOUBLE, (void*)(*arg2));
      if (!array) SWIG_fail;
      resultobj = SWIG_Python_AppendOutput(resultobj,array);
    }
    return resultobj;
  fail:
    return NULL;
  }
  %}

%typemap(in) (const char *name)
{
  $1 = PyString_AsString($input);
}
%init %{
  import_array();
  %}

%include "cow.h"
%native(getarray3) PyObject *_wrap_getarray3(PyObject *SWIGUNUSEDPARM(self), PyObject *args);

%callback("%(upper)s");
void test_trans(double *result, double **args, int **strides, void *udata);
                   %nocallback;
