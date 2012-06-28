


%module cow

%typemap(in) (const char *name)
{
  $1 = PyString_AsString($input);
}

%inline %{
#include "cow.h"
%}
%include "cow.h"
