
from cowpy cimport *
from libc.stdlib cimport malloc, free
import numpy as np

cdef extern from "srhdpack.h":
    int SRHDPACK_VELOCITY_GAMMA          = -42
    int SRHDPACK_VELOCITY_BETA           = -43
    int SRHDPACK_VELOCITY_GAMMABETA      = -44
    int SRHDPACK_VELOCITY_DUMUDXMU       = -45
    int SRHDPACK_VELOCITY_DUMUDUMU       = -46
    int SRHDPACK_SEPARATION_LAB          = -47
    int SRHDPACK_SEPARATION_PROPER       = -48
    int SRHDPACK_PROJECTION_NONE         = -49
    int SRHDPACK_PROJECTION_TRANSVERSE   = -50
    int SRHDPACK_PROJECTION_LONGITUDINAL = -51

    struct srhdpack_samplemode:
        double exponent # exponent value: p
        int velmode     # SRHDPACK_VELOCITY
        int sepmode     # SRHDPACK_SEPARATION
        int projmode    # SRHDPACK_PROJECTION
    
    void srhdpack_shelevequescaling(cow_dfield *vel,
                                    cow_histogram *hist,
                                    int velmode,
                                    int sepmode,
                                    int projmode,
                                    int nbatch,
                                    int nperbatch,
                                    int seed,
                                    double exponent)

    void srhdpack_collectpairs(cow_dfield *vel,
                               srhdpack_samplemode *modes,
                               int num_modes,
                               int num_pairs,
                               int num_samps,
                               double *samploc,
                               double *outbufx,
                               double *outbufy)


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



def inverse_dict(d):
    return dict((v,k) for k, v in d.iteritems())

_velocity   = {"gamma"        : SRHDPACK_VELOCITY_GAMMA,
               "beta"         : SRHDPACK_VELOCITY_BETA,
               "gammabeta"    : SRHDPACK_VELOCITY_GAMMABETA,
               "dumudxmu"     : SRHDPACK_VELOCITY_DUMUDXMU,
               "dumudumu"     : SRHDPACK_VELOCITY_DUMUDUMU}
_separation = {"lab"          : SRHDPACK_SEPARATION_LAB,
               "proper"       : SRHDPACK_SEPARATION_PROPER}
_projection = {"none"         : SRHDPACK_PROJECTION_NONE,
               "transverse"   : SRHDPACK_PROJECTION_TRANSVERSE,
               "longitudinal" : SRHDPACK_PROJECTION_LONGITUDINAL}
_velocity_i = inverse_dict(_velocity)
_separation_i = inverse_dict(_separation)
_projection_i = inverse_dict(_projection)


cdef class PairwiseStructureFunction(object):
    cdef double x0, x1
    cdef srhdpack_samplemode _c

    def __init__(self, **kwargs):
        self._c.exponent = 1.0
        self._c.velmode = SRHDPACK_VELOCITY_GAMMABETA
        self._c.sepmode = SRHDPACK_SEPARATION_PROPER
        self._c.projmode = SRHDPACK_PROJECTION_LONGITUDINAL
        for k, v in kwargs.iteritems():
            setattr(self, k, v)

    property exponent:
        def __get__(self):
            return self._c.exponent
        def __set__(self, exponent):
            self._c.exponent = exponent

    property velocity:
        def __get__(self):
            return _velocity_i[self._c.velmode]
        def __set__(self, mode):
            self._c.velmode = _velocity[mode]

    property separation:
        def __get__(self):
            return _separation_i[self._c.sepmode]
        def __set__(self, mode):
            self._c.sepmode = _separation[mode]

    property projection:
        def __get__(self):
            return _projection_i[self._c.projmode]
        def __set__(self, mode):
            self._c.projmode = _projection[mode]


def structure_functions(DataField vel, struc_fns, samples=100, reuse=1,
                        maxmem=10, x0=1e-3, x1=1.0, spacing="linear", bins=128):
    """
    Description
    -----------

    Processes batches of pair samples from the distributed vector field
    `vel` and returns the pairs processed by each object in `struc_fns`.


    Parameters
    ----------

    samples ... Number of total field locations to be sampled across all procs.

    reuse ..... Average number of times (per processing mode) each sample will
                be used during pair processsing. The default value of 1
                corresponds to taking `samples` / 2 pairs.

    maxmem .... Memory, in MB, allowed to be used. Sampling will be broken into
                batches if maxmem is less than the approximate required to run
                the sampling.

    x0 ........ Lower bound on histogram.

    x1 ........ Upper bound on histogram.

    spacing ... Either 'linear' or 'log'.

    bins ...... Number of histogram bins.

    """

    assert vel is not None

    # 
    # Each proc will take its own share of the samples and pairs
    # 
    cdef long num_samps = samples / vel.domain.cart_size
    cdef long num_pairs = num_samps / 2 * reuse
    cdef long num_modes = len(struc_fns)

    #
    # Approximate the memory usage on each process
    #
    cdef long bytes = (num_samps*3 + num_pairs*num_modes*2) * sizeof(double)
    cdef long batches = bytes / (maxmem * 1024 * 1024) + 1
    num_samps /= batches
    num_pairs /= batches

    print "processing request for %d total pairs from %d total samples" % (
        samples * reuse / 2, samples)
    print "total memory requirement per process is %d MB and maxmem is %d MB, "\
        "running %d batches each with %d samples" % (bytes / (1024 * 1024),
                                                     maxmem, batches, num_samps)
    print "seeding Numpy's random number generator with cartesian rank"
    np.random.seed(vel.domain.cart_rank)

    cdef PairwiseStructureFunction S
    cdef long m
    cdef srhdpack_samplemode *modes = <srhdpack_samplemode*> malloc(
        num_modes * sizeof(srhdpack_samplemode))

    cdef np.ndarray[np.double_t,ndim=2] samploc = np.zeros([num_samps, 3])
    cdef np.ndarray[np.double_t,ndim=2] outbufx = np.zeros([num_pairs, num_modes])
    cdef np.ndarray[np.double_t,ndim=2] outbufy = np.zeros([num_pairs, num_modes])

    hists = [Histogram1d(x0, x1, bins=bins, spacing=spacing, domain=vel.domain,
                         binmode="average") for s in struc_fns]

    for m, S in enumerate(struc_fns):
        modes[m] = S._c
    for m in range(batches):
        print "running batch %d/%d" % (m+1, batches)
        samploc[...] = np.random.rand(num_samps, 3)
        srhdpack_collectpairs(vel._c, modes, num_modes, num_pairs, num_samps,
                              <double*>samploc.data,
                              <double*>outbufx.data,
                              <double*>outbufy.data)
        for n, hist in enumerate(hists):
            for xval, yval in zip(outbufx[:,n], outbufy[:,n]):
                hist.add_sample(xval, yval)
    free(modes)
    for hist in hists:
        hist.seal()
    return hists
