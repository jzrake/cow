


%module cow

%typemap(in) (const char *name)
{
  $1 = PyString_AsString($input);
}

%inline %{
  extern "C"{
#include "cow.h"
#include "histogram.hpp"
  }
%}
%include "cow.h"
%include "histogram.hpp"
