#!/usr/bin/env python

from distutils.core import setup, Extension
import numpy as np

# Helpful info on linker environment:
# http://bit.ly/Nsvato

config = {
    'COW_HDF5': 0,
    'COW_HDF5_MPI': 0,
    'COW_FFTW': 0,
    'COW_MPI': 0,
    'include_dirs': [ ],
    'library_dirs': [ ],
    'libraries': [ ], # e.g. ['hdf5', 'z', 'fftw']
    'extra_compile_args': [ ],
    'extra_link_args': [ ],
    'NPY_INC': np.get_include() }

try:
    import cow_config
    config.update(cow_config.config)
    print "Using system config"
except:
    print "No system config, using default settings"
config['include_dirs'] += [np.get_include()]

csource = ['cow.c', 'io.c', 'hist.c', 'samp.c', 'fft.c', 'fft_3d.c',
           'remap_3d.c', 'pack_3d.c']

cow_module = Extension('_cow',
        extra_compile_args = ['-std=c99'] + config['extra_compile_args'],
        extra_link_args = config['extra_link_args'],
        define_macros = [a for a in config.items() if a[0].startswith('COW')],
        include_dirs = ["../src"] + config['include_dirs'],
        library_dirs = config['library_dirs'],
        libraries = config['libraries'],
        sources = ['cow.i'] + ["../src/" + c for c in csource])

setup(name        = 'cow',
      version     = '0.4',
      author      = "Jonathan Zrake",
      description = """C.O.W.""",
      ext_modules = [cow_module],
      py_modules  = ["cowpy"])
