
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

def relative_lorentz_pairs(DataField vel, nsamples, bins=36, x0=1e-3, x1=1.0,
                           nperbatch=10000, seed=True, spacing="linear"):
    histpro = Histogram1d(x0, x1, bins=bins, spacing=spacing, commit=False)
    histlab = Histogram1d(x0, x1, bins=bins, spacing=spacing, commit=False)
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


def sheleveque_scaling(DataField vel, nsamples, maxp=10, bins=36, x0=1e-3, x1=1.0,
                       nperbatch=10000, seed=True, spacing="linear"):
    hists = { }
    for p in range(1,maxp+1):
        hist = Histogram1d(x0, x1, bins=bins, spacing=spacing, commit=False)
        srhdpack_shelevequescaling(vel._c, hist._c,
                                   SRHDPACK_VELOCITY_GAMMABETA,
                                   SRHDPACK_SEPARATION_PROPER,
                                   SRHDPACK_PROJECTION_LONGITUDINAL,
                                   int(nsamples/nperbatch), nperbatch, seed, 0.5*p)
        hists[p] = hist
    return hists
