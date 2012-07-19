
from capi.ccow import srhdpack_relativelorentzpairs
from cowpy import Histogram1d

def relative_lorentz_pairs(vel, nsamples, nperbatch=10000, seed=True):
    histpro = Histogram1d(0.1, 1.5, spacing="linear", commit=False)
    histlab = Histogram1d(0.1, 1.5, spacing="linear", commit=False)
    srhdpack_relativelorentzpairs(vel._c, histpro._c, histlab._c,
                                  int(nsamples/nperbatch), nperbatch, seed)
    return histpro, histlab
