
from cowpy cimport *

cdef extern from "srhdpack.h":
    cdef int SRHDPACK_VELOCITY_GAMMA          = -42
    cdef int SRHDPACK_VELOCITY_BETA           = -43
    cdef int SRHDPACK_VELOCITY_GAMMABETA      = -44
    cdef int SRHDPACK_SEPARATION_LAB          = -45
    cdef int SRHDPACK_SEPARATION_PROPER       = -46
    cdef int SRHDPACK_PROJECTION_NONE         = -47
    cdef int SRHDPACK_PROJECTION_TRANSVERSE   = -49
    cdef int SRHDPACK_PROJECTION_LONGITUDINAL = -48

    void srhdpack_shelevequescaling(cow_dfield *vel,
                                    cow_histogram *hist,
                                    int velmode,
                                    int sepmode,
                                    int projmode,
                                    int nbatch,
                                    int nperbatch,
                                    int seed,
                                    double exponent)

def relative_lorentz_pairs(DataField vel, nsamples, bins=36, nperbatch=10000,
                           seed=True):
    histpro = Histogram1d(0.0, 1.0, bins=bins, spacing="linear", commit=False)
    histlab = Histogram1d(0.0, 1.0, bins=bins, spacing="linear", commit=False)
    srhdpack_shelevequescaling(vel._c, histpro._c,
                               SRHDPACK_VELOCITY_GAMMA,
                               SRHDPACK_SEPARATION_PROPER,
                               SRHDPACK_PROJECTION_NONE,
                               int(nsamples/nperbatch), nperbatch, seed, 1.0)
    srhdpack_shelevequescaling(vel._c, histlab._c,
                               SRHDPACK_VELOCITY_GAMMA,
                               SRHDPACK_SEPARATION_LAB,
                               SRHDPACK_PROJECTION_NONE,
                               int(nsamples/nperbatch), nperbatch, seed, 1.0)
    return histpro, histlab


def sheleveque_scaling(DataField vel, nsamples, maxp=10, bins=36,
                       nperbatch=10000, seed=True):
    hists = { }
    for p in range(1,maxp+1):
        hist = Histogram1d(0.0, 1.0, bins=bins, spacing="linear", commit=False)
        srhdpack_shelevequescaling(vel._c, hist._c,
                                   SRHDPACK_VELOCITY_GAMMA,
                                   SRHDPACK_SEPARATION_PROPER,
                                   SRHDPACK_PROJECTION_NONE,
                                   int(nsamples/nperbatch), nperbatch, seed, 0.5*p)
        hists[p] = hist
    return hists
