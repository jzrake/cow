
from cowpy cimport *
from cowpy import DataField, Histogram1d

cdef extern from "srhdpack.h":
    void srhdpack_relativelorentzpairs(cow_dfield *vel,
                                       cow_histogram *histpro,
                                       cow_histogram *histlab,
                                       int nbatch,
                                       int nperbatch,
                                       int seed)

def relative_lorentz_pairs(DataField vel, nsamples, bins=36, nperbatch=10000, seed=True):
    histpro = Histogram1d(0.0, 2.5, bins=bins, spacing="linear", commit=False)
    histlab = Histogram1d(0.0, 2.5, bins=bins, spacing="linear", commit=False)
    srhdpack_relativelorentzpairs(vel._c, histpro._c, histlab._c,
                                  int(nsamples/nperbatch), nperbatch, seed)
    return histpro, histlab
