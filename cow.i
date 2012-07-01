


%module cow

%typemap(in) (const char *name)
{
  $1 = PyString_AsString($input);
}


%inline %{
#include "cow.h"
int test_trans(double *result, double **args, int **strides, void *udata);
%}

%include "cow.h"
%constant int test_trans(double *result, double **args, int **strides, void *udata);


