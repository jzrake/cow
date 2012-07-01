


%module cow

%typemap(in) (const char *name)
{
  $1 = PyString_AsString($input);
}

%inline %{
  extern "C"{
#include "cow.h"
  }
%}
%include "cow.h"

