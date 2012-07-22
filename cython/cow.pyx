
cimport numpy as np

cdef extern from "cow.h":
    cdef enum:
        COW_NOREOPEN_STDOUT = (2<<0)
        COW_NOREOPEN_STDOUT = (2<<0)
        COW_DISABLE_MPI     = (2<<1)
        COW_ALL_DIMS             = -41
        COW_HIST_SPACING_LINEAR  = -42
        COW_HIST_SPACING_LOG     = -43
        COW_HIST_BINMODE_COUNTS  = -44 # traditional histogram
        COW_HIST_BINMODE_DENSITY = -45 # divides by bin width
        COW_HIST_BINMODE_AVERAGE = -46
        COW_PROJECT_OUT_DIV      = -47 # used for Helmholtz decomposition
        COW_PROJECT_OUT_CURL     = -48
        COW_SAMPLE_NEAREST       = -49 # sample the nearest zone center
        COW_SAMPLE_LINEAR        = -50 # use (uni/bi/tri) linear interp
        COW_SAMPLE_ERROR_OUT     = -51 # out-of-bounds sample request
        COW_SAMPLE_ERROR_WRONGD  = -52 # wrong number of dims on sample coords

    cdef struct cow_domain
    cdef struct cow_dfield
    cdef struct cow_histogram
    ctypedef void (*cow_transform)(double *result, double **args, int **strides,
                                   void *udata)

    void cow_init(int argc, char **argv, int modes)
    void cow_finalize()
    int cow_mpirunning()

    cow_domain *cow_domain_new()
    void cow_domain_commit(cow_domain *d)
    void cow_domain_del(cow_domain *d)
    void cow_domain_setsize(cow_domain *d, int dim, int size)
    void cow_domain_setndim(cow_domain *d, int ndim)
    void cow_domain_setguard(cow_domain *d, int guard)
    void cow_domain_setprocsizes(cow_domain *d, int dim, int size)
    void cow_domain_setcollective(cow_domain *d, int mode)
    void cow_domain_setchunk(cow_domain *d, int mode)
    void cow_domain_setalign(cow_domain *d, int alignthreshold, int diskblocksize)
    void cow_domain_readsize(cow_domain *d, char *fname, char *dname)
    int cow_domain_getndim(cow_domain *d)
    int cow_domain_getguard(cow_domain *d)
    int cow_domain_getnumlocalzonesincguard(cow_domain *d, int dim)
    int cow_domain_getnumlocalzonesinterior(cow_domain *d, int dim)
    int cow_domain_getnumglobalzones(cow_domain *d, int dim)
    int cow_domain_getglobalstartindex(cow_domain *d, int dim)
    int cow_domain_getgridspacing(cow_domain *d, int dim)
    int cow_domain_getcartrank(cow_domain *d)
    int cow_domain_getcartsize(cow_domain *d)
    int cow_domain_subgridatposition(cow_domain *d, double x, double y, double z)
    int cow_domain_indexatposition(cow_domain *d, int dim, double x)
    double cow_domain_positionatindex(cow_domain *d, int dim, int index)
    void cow_domain_barrier(cow_domain *d)

    cow_dfield *cow_dfield_new(cow_domain *domain, char *name)
    cow_dfield *cow_dfield_dup(cow_dfield *f)
    void cow_dfield_commit(cow_dfield *f)
    void cow_dfield_del(cow_dfield *f)
    void cow_dfield_addmember(cow_dfield *f, char *name)
    void cow_dfield_setname(cow_dfield *f, char *name)
    void cow_dfield_extract(cow_dfield *f, int *I0, int *I1, void *out)
    void cow_dfield_replace(cow_dfield *f, int *I0, int *I1, void *out)
    void cow_dfield_loop(cow_dfield *f, cow_transform op, void *udata)
    void cow_dfield_settransform(cow_dfield *f, cow_transform op)
    void cow_dfield_clearargs(cow_dfield *f)
    void cow_dfield_pusharg(cow_dfield *f, cow_dfield *arg)
    void cow_dfield_setuserdata(cow_dfield *f, void *userdata)
    void cow_dfield_setiparam(cow_dfield *f, int p)
    void cow_dfield_setfparam(cow_dfield *f, double p)
    void cow_dfield_transformexecute(cow_dfield *f)
    char *cow_dfield_iteratemembers(cow_dfield *f)
    char *cow_dfield_nextmember(cow_dfield *f)
    char *cow_dfield_getname(cow_dfield *f)
    cow_domain *cow_dfield_getdomain(cow_dfield *f)
    int cow_dfield_getstride(cow_dfield *f, int dim)
    int cow_dfield_getnmembers(cow_dfield *f)
    size_t cow_dfield_getdatabytes(cow_dfield *f)
    void cow_dfield_setbuffer(cow_dfield *f, void *buffer)
    void cow_dfield_sampleglobalind(cow_dfield *f, int i, int j, int k, double **x, int *n0)
    int cow_dfield_setsamplecoords(cow_dfield *f, double *x, int n0, int n1)
    void cow_dfield_getsamplecoords(cow_dfield *f, double **x, int *n0, int *n1)
    void cow_dfield_getsampleresult(cow_dfield *f, double **x, int *n0, int *n1)
    void cow_dfield_setsamplemode(cow_dfield *f, int mode)
    void cow_dfield_sampleexecute(cow_dfield *f)
    int cow_dfield_getownsdata(cow_dfield *f)
    void *cow_dfield_getbuffer(cow_dfield *f)
    void cow_dfield_syncguard(cow_dfield *f)
    void cow_dfield_reduce(cow_dfield *f, double x[3])
    void cow_dfield_write(cow_dfield *f, char *fname)
    void cow_dfield_read(cow_dfield *f, char *fname)

    cow_histogram *cow_histogram_new()
    void cow_histogram_commit(cow_histogram *h)
    void cow_histogram_del(cow_histogram *h)
    void cow_histogram_setbinmode(cow_histogram *h, int binmode)
    void cow_histogram_setspacing(cow_histogram *h, int spacing)
    void cow_histogram_setnbins(cow_histogram *h, int dim, int nbinsx)
    void cow_histogram_setlower(cow_histogram *h, int dim, double v0)
    void cow_histogram_setupper(cow_histogram *h, int dim, double v1)
    void cow_histogram_setfullname(cow_histogram *h, char *fullname)
    void cow_histogram_setnickname(cow_histogram *h, char *nickname)
    void cow_histogram_setdomaincomm(cow_histogram *h, cow_domain *d)
    void cow_histogram_addsample1(cow_histogram *h, double x, double w)
    void cow_histogram_addsample2(cow_histogram *h, double x, double y, double w)
    void cow_histogram_dumpascii(cow_histogram *h, char *fn)
    void cow_histogram_dumphdf5(cow_histogram *h, char *fn, char *dn)
    void cow_histogram_seal(cow_histogram *h)
    int cow_histogram_getsealed(cow_histogram *h)
    long cow_histogram_gettotalcounts(cow_histogram *h)
    void cow_histogram_populate(cow_histogram *h, cow_dfield *f, cow_transform op)
    void cow_histogram_getbinlocx(cow_histogram *h, double **x, int *n0)
    void cow_histogram_getbinlocy(cow_histogram *h, double **x, int *n0)
    void cow_histogram_getbinval1(cow_histogram *h, double **x, int *n0)
    void cow_histogram_getbinval2(cow_histogram *h, double **x, int *n0, int *n1)
    double cow_histogram_getbinval(cow_histogram *h, int i, int j)
    char *cow_histogram_getname(cow_histogram *h)

    void cow_fft_pspecvecfield(cow_dfield *f, cow_histogram *h)
    void cow_fft_helmholtzdecomp(cow_dfield *f, int mode)

    void cow_trans_divcorner(double *result, double **args, int **s, void *u)
    void cow_trans_div5(double *result, double **args, int **s, void *u)
    void cow_trans_rot5(double *result, double **args, int **s, void *u)
    void cow_trans_component(double *result, double **args, int **s, void *u)
    void cow_trans_magnitude(double *result, double **args, int **s, void *u)
    void cow_trans_cross(double *result, double **args, int **s, void *u)
    void cow_trans_dot3(double *result, double **args, int **s, void *u)


import os
import sys
import atexit
import warnings
import numpy as np
import h5py
try:
    import h5py
except ImportError:
    warnings.warn("h5py was not detected, you might be missing some functions")



def _getie(s):
    return int(os.getenv(s, 0))
def exitfunc():
    cow_finalize()

cdef _init_cow():
    modes = 0
    _hdf5_collective = _getie("COW_HDF5_COLLECTIVE")
    modes |= (COW_NOREOPEN_STDOUT if _getie("COW_NOREOPEN_STDOUT") else 0)
    modes |= (COW_DISABLE_MPI if _getie("DISABLE_MPI") else 0)
    modes |= (COW_DISABLE_MPI if '-s' in sys.argv else 0)
    cdef int argc = 0
    cdef char *argv[1]
    cow_init(argc, argv, modes)
    atexit.register(exitfunc)

_init_cow()
_hdf5_collective = 0


cdef class DistributedDomain(object):
    cdef cow_domain *_c
    def __cinit__(self, G_ntot, guard=0, *args, **kwargs):
        cdef int KILOBYTES = 1 << 10
        cdef int MEGABYTES = 1 << 20
        print "building domain", G_ntot
        nd = len(G_ntot)
        assert nd <= 3
        self._c = cow_domain_new()
        for n, ni in enumerate(G_ntot):
            cow_domain_setsize(self._c, n, ni)
        cow_domain_setndim(self._c, nd)
        cow_domain_setguard(self._c, guard)
        cow_domain_commit(self._c)
        cow_domain_setchunk(self._c, 1) # set up the IO scheme after commit
        cow_domain_setcollective(self._c, _hdf5_collective)
        cow_domain_setalign(self._c, 4*KILOBYTES, 4*MEGABYTES)

    def __dealloc__(self):
        if self._c: # in case __cinit__ raised something, don't clean this up
            cow_domain_del(self._c)

    @property
    def cart_rank(self):
        return cow_domain_getcartrank(self._c)

    @property
    def cart_size(self):
        return cow_domain_getcartsize(self._c)

    @property
    def global_start(self):
        return tuple([cow_domain_getglobalstartindex(self._c, n)
                      for n in range(self.ndim)])

    @property
    def global_shape(self):
        return tuple([cow_domain_getnumglobalzones(self._c, n)
                      for n in range(self.ndim)])
    @property
    def ndim(self):
        return cow_domain_getndim(self._c)

    @property
    def guard(self):
        return cow_domain_getguard(self._c)

    def coordinate(self, ind):
        assert len(ind) == self.ndim
        return tuple([cow_domain_positionatindex(self._c, n, i)
                      for n, i in enumerate(ind)])

    def barrier(self):
        cow_domain_barrier(self._c)

    def sequential(self, func, args=()):
        for i in range(self.cart_size):
            if self.cart_rank == i:
                func(*args)
            self.barrier()


cdef class DataField(object):
    cdef cow_dfield *_c
    cdef np.ndarray _buf
    cdef DistributedDomain _domain
    def __cinit__(self, DistributedDomain domain, members, name="datafield", *args, **kwargs):
        self._c = cow_dfield_new(domain._c, name)
        self._domain = domain
        for m in members:
            cow_dfield_addmember(self._c, m)
        nd = domain.ndim
        dims = [ ]
        for i in range(nd):
            dims.append(cow_domain_getnumlocalzonesincguard(domain._c, i))
        dims.append(len(members))

        cdef np.ndarray[np.double_t,ndim=2] _buf1
        cdef np.ndarray[np.double_t,ndim=3] _buf2
        cdef np.ndarray[np.double_t,ndim=4] _buf3

        if nd == 1:
            _buf1 = np.zeros(dims)
            self._buf = _buf1
            cow_dfield_setbuffer(self._c, <double*>_buf1.data)
        elif nd == 2:
            _buf2 = np.zeros(dims)
            self._buf = _buf2
            cow_dfield_setbuffer(self._c, <double*>_buf2.data)
        elif nd == 3:
            _buf3 = np.zeros(dims)
            self._buf = _buf3
            cow_dfield_setbuffer(self._c, <double*>_buf3.data)

        cow_dfield_commit(self._c)

        #self._lookup = dict([(m, n) for n,m in enumerate(members)])

    def __dealloc__(self):
        if self._c: # in case __cinit__ raised something, don't clean this up
            cow_dfield_del(self._c)

    @property
    def domain(self):
        return self._domain

    property name:
        def __get__(self):
            return cow_dfield_getname(self._c)
        def __set__(self, name):
            cow_dfield_setname(self._c, name)

    property members:
        def __get__(self):
            mem = [cow_dfield_iteratemembers(self._c)]
            for n in range(cow_dfield_getnmembers(self._c) - 1):
                mem.append(cow_dfield_nextmember(self._c))
            return mem

    property value:
        def __get__(self):
            return self._buf
        def __set__(self, val):
            self._buf[...] = val

    property interior:
        """
        Opertates on a view of the buffer object which refers only to the
        interior portion, excluding the guard zones.
        """
        def __get__(self):
            ng = self.domain.guard
            nd = self.domain.ndim
            if ng == 0:
                return self._buf
            elif nd == 1:
                return self._buf[ng:-ng, :]
            elif nd == 2:
                return self._buf[ng:-ng, ng:-ng, :]
            elif nd == 3:
                return self._buf[ng:-ng, ng:-ng, ng:-ng, :]
        def __set__(self, val):
            ng = self.domain.guard
            nd = self.domain.ndim
            if ng == 0:
                self._buf[:] = val
            elif nd == 1:
                self._buf[ng:-ng, :] = val
            elif nd == 2:
                self._buf[ng:-ng, ng:-ng, :] = val
            elif nd == 3:
                self._buf[ng:-ng, ng:-ng, ng:-ng, :] = val

    cdef sync_guard(self):
        cow_dfield_syncguard(self._c)

    cdef dump(self, fname):
        cow_dfield_write(self._c, fname)

    cdef read(self, fname):
        cow_dfield_read(self._c, fname)

    cdef _apply_transform(self, args, cow_transform op, void *userdata=NULL):
        """
        Private method, removes redundant code in the two functions that follow.
        """
        cow_dfield_clearargs(self._c)
        for arg in args:
            cow_dfield_pusharg(self._c, <cow_dfield*>arg._c)
        cow_dfield_settransform(self._c, op)
        cow_dfield_setuserdata(self._c, userdata)
        cow_dfield_transformexecute(self._c)
        return self

    cdef reduce_component(self, member):
        """
        Returns the min, max, and sum of the data member `member`, which may be
        int or str.
        """
        if type(member) is str:
            member = self._lookup[member]
        else:
            assert member < len(self._members)
        cow_dfield_clearargs(self._c)
        cow_dfield_settransform(self._c, cow_trans_component)
        cow_dfield_setuserdata(self._c, self._c)
        cow_dfield_setiparam(self._c, member)
        cdef np.ndarray[np.double_t,ndim=1] res = np.zeros(3)
        cow_dfield_reduce(self._c, <double*>res.data)
        return res
    
    cdef reduce_magnitude(self):
        """
        Returns the min, max, and sum of the data fields's vector magnitude.
        """
        cow_dfield_clearargs(self._c)
        cow_dfield_settransform(self._c, cow_trans_magnitude)
        cow_dfield_setuserdata(self._c, self._c)
        cdef np.ndarray[np.double_t,ndim=1] res = np.zeros(3)
        cow_dfield_reduce(self._c, <double*>res.data)
        return res

    def __getitem__(self, key):
        if type(key) is int:
            return self.value[..., key]
        else:
            return self.value[..., self._lookup[key]]

    def __setitem__(self, key, val):
        if type(key) is int:
            self.value[..., key] = val 
        else:
            self.value[..., self._lookup[key]] = val

    def __repr__(self):
        props = ["%s" % type(self),
                 "name: %s" % self.name,
                 "members: %s" % str(self.members),
                 "local shape: %s" % str(self.value[..., 0].shape),
                 "global shape: %s" % str(self.domain.global_shape),
                 "global start: %s" % str(self.domain.global_start),
                 "padding: %s" % self.domain.guard]
        return "{" + "\n\t".join(props) + "}"
