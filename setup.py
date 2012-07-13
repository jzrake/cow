#!/usr/bin/env python

from distutils.core import setup, Extension
import numpy as np

config = {
    'COW_HDF5': 0,
    'COW_HDF5_MPI': 0,
    'COW_FFTW': 0,
    'COW_MPI': 0,
    'include_dirs': [ ],
    'library_dirs': [ ],
    'libraries': [ ], # e.g. ['hdf5', 'z', 'fftw']
    'NPY_INC': np.get_include() }

try:
    import cow_config
    config.update(cow_config.config)
    print "Using system config"
except:
    print "No system config, using default settings"
config['include_dirs'] += [np.get_include()]

cow_module = Extension \
('_cow',
 extra_compile_args=['-std=c99'],
 define_macros = [a for a in config.items() if a[0].startswith('COW')],
 include_dirs = config['include_dirs'],
 library_dirs = config['library_dirs'],
 libraries = config['libraries'],
 sources = ['cow.i',
            'cow.c',
            'io.c',
            'hist.c',
            'samp.c',
            'fft.c',
            'fft_3d.c',
            'remap_3d.c',
            'pack_3d.c',
            'factor.c'])
setup(name        = 'cow',
      version     = '0.4',
      author      = "Jonathan Zrake",
      description = """C.O.W.""",
      ext_modules = [cow_module],
      py_modules  = ["cow"])
