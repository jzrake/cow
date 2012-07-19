/* -*- C -*-  (not really, but good for syntax highlighting) */

%module csrhdpack
%{

#define SWIG_FILE_WITH_INIT

  // Swig source is littered with 'no prototype' warnings on icc:
#ifdef __INTEL_COMPILER
#pragma warning disable 1418
#endif

  //#include <numpy/arrayobject.h>
#include "../src/cow.h"

void srhdpack_relativelorentzpairs(cow_dfield *vel,
				   cow_histogram *histpro,
				   cow_histogram *histlab,
				   int nbatch,
				   int nperbatch,
				   int seed);
  %}

void srhdpack_relativelorentzpairs(cow_dfield *vel,
				   cow_histogram *histpro,
				   cow_histogram *histlab,
				   int nbatch,
				   int nperbatch,
				   int seed);
