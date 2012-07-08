#!/usr/bin/env python

from distutils.core import setup, Extension


cow_module = Extension \
('_cow',
 extra_compile_args=['-std=c99'],
 define_macros = [('COW_MPI', '1'),
                  ('COW_HDF5', '1'),
                  ('COW_FFTW', '1'),
                  ('FFT_FFTW', '1'),
                  ('COW_HDF5_MPI', '1')],
 include_dirs = ['/Library/Science/hdf5-1.8.9-par/include',
                 '/Library/Science/fftw-2.1.5/include',
                 '/Library/Science/mpich2/include',
                 '/Library/Frameworks/Python.framework/Versions/7.1/lib/python2.7/site-packages/numpy/core/include'],
 library_dirs = ['/Library/Science/hdf5-1.8.9-par/lib',
                 '/Library/Science/fftw-2.1.5/lib',
                 '/Library/Science/mpich2/lib'],
 libraries = ['hdf5', 'z', 'fftw', 'mpich', 'pmpich'],
 sources=['cow.i',
          'cow.c',
          'io.c',
          'hist.c',
          'fft.c',
          'fft_3d.c',
          'remap_3d.c',
          'pack_3d.c',
          'factor.c'])
setup(name        = 'cow',
      version     = '0.1',
      author      = "Jonathan Zrake",
      description = """C.O.W.""",
      ext_modules = [cow_module],
      py_modules  = ["cow"])
