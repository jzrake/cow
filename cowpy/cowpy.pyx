
cimport numpy as np
cimport cowpy

import os
import sys
import atexit
import warnings
import numpy as np
try:
    import h5py
except ImportError:
    warnings.warn("h5py was not detected, you might be missing some functions")



def _getie(s, dflt=0):
    return int(os.getenv(s, dflt))
def exitfunc():
    cow_finalize()


cdef _init_cow():
    modes = 0
    _runtime_cfg['hdf5_collective'] = _getie("COW_HDF5_COLLECTIVE")
    _runtime_cfg['hdf5_chunk'] = _getie("COW_HDF5_CHUNK", dflt=1)
    modes |= (COW_NOREOPEN_STDOUT if _getie("COW_NOREOPEN_STDOUT") else 0)
    modes |= (COW_DISABLE_MPI if _getie("COW_DISABLE_MPI") else 0)
    cdef int argc = 0
    cdef char *argv[1]
    cow_init(argc, argv, modes)
    atexit.register(exitfunc)

_runtime_cfg = {'hdf5_collective': 0}
_init_cow()


cdef class DistributedDomain(object):
    def __cinit__(self):
        self._c = cow_domain_new()

    def __init__(self, G_ntot, guard=0, *args, **kwargs):
        cdef int KILOBYTES = 1 << 10
        cdef int MEGABYTES = 1 << 20
        print "building domain", G_ntot
        nd = len(G_ntot)
        if nd > 3:
            raise ValueError("domain dims must be no larger 3")
        for n, ni in enumerate(G_ntot):
            cow_domain_setsize(self._c, n, ni)
        cow_domain_setndim(self._c, nd)
        cow_domain_setguard(self._c, guard)
        cow_domain_commit(self._c)
        # set up the IO scheme after commit
        cow_domain_setchunk(self._c, _runtime_cfg['hdf5_chunk'])
        cow_domain_setcollective(self._c, _runtime_cfg['hdf5_collective'])
        cow_domain_setalign(self._c, 4*KILOBYTES, 4*MEGABYTES)

    def __dealloc__(self):
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
    def __cinit__(self):
        self._c = cow_dfield_new()

    def __init__(self, DistributedDomain domain, members=[], name="datafield",
                 *args, **kwargs):
        if domain is None or len(members) == 0:
            raise ValueError("bad argument list")

        cow_dfield_setdomain(self._c, domain._c)
        cow_dfield_setname(self._c, name)
        self._domain = domain

        for m in members:
            cow_dfield_addmember(self._c, m)
        nd = domain.ndim
        dims = [ ]
        for i in range(nd):
            dims.append(cow_domain_getnumlocalzonesincguard(domain._c, i))
        dims.append(len(members))

        if nd == 1:
            self._buf = np.zeros(dims)
            self._flg = np.zeros(dims[0:1], dtype=np.int32)
            cow_dfield_setdatabuffer(self._c, <double*>self._buf.data)
            cow_dfield_setflagbuffer(self._c, <int*>self._flg.data)
        elif nd == 2:
            self._buf = np.zeros(dims)
            self._flg = np.zeros(dims[0:2], dtype=np.int32)
            cow_dfield_setdatabuffer(self._c, <double*>self._buf.data)
            cow_dfield_setflagbuffer(self._c, <int*>self._flg.data)
        elif nd == 3:
            self._buf = np.zeros(dims)
            self._flg = np.zeros(dims[0:3], dtype=np.int32)
            cow_dfield_setdatabuffer(self._c, <double*>self._buf.data)
            cow_dfield_setflagbuffer(self._c, <int*>self._flg.data)
        cow_dfield_commit(self._c)

    def __dealloc__(self):
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

    property flags:
        def __get__(self):
            return self._flg
        
        def __set__(self, val):
            self._flg[...] = val

    def sync_guard(self):
        cow_dfield_syncguard(self._c)

    def dump(self, fname):
        cow_dfield_write(self._c, fname)

    def read(self, fname):
        cow_dfield_read(self._c, fname)

    def sample_global(self, pnts):
        """
        Samples the distributed domain at the provided list of physical
        coordinates, `pnts` which has the shape (N x nd). Returns a permutation
        of those coordinates, and their associated field values. Values are
        multi-linearly interpolated from cell centers. This function must be
        called by all ranks which own this array, although the size of `pnts`
        may vary between ranks, and is allowed to be zero.

        NOTE: Typically this function is used for the calculation of structure
        functions over various fields. Due to the choice of algorithm for
        distributing the data among parallel processes, the list of output
        points is sorted according to the subgrid containing the target
        coordinate, and thus contains spatial correlations in the ordering of
        the output points. This means that a random sampling of pairs (or
        triplets, etc.) requires the user to randomly permute the output vector,
        or to index it with randomly selected indices.
        """
        cdef np.ndarray[np.double_t,ndim=2] points = np.array(pnts)
        if points.shape[1] != self.domain.ndim:
            raise ValueError("dimension of sample coordinates must match those "
                             "of the domain")
        if self.domain.guard <= 0:
            raise ValueError("sampling requires at least one guard zone")
        cdef np.ndarray[np.double_t,ndim=2] p3d = np.zeros([points.shape[0], 3])
        for n in range(self.domain.ndim):
            p3d[:,n] = points[:,n]
        err = cow_dfield_setsamplecoords(self._c, <double*>p3d.data,
                                         p3d.shape[0], 3)
        if err != 0:
            raise RuntimeError("off-bounds sample coordinates")
        cow_dfield_setsamplemode(self._c, COW_SAMPLE_LINEAR)
        cow_dfield_sampleexecute(self._c)
        cdef double *x
        cdef double *P
        cdef int nd = self.domain.ndim
        cdef int nq
        cdef int Nx
        cdef int NP
        cow_dfield_getsamplecoords(self._c, &x, &Nx, NULL)
        cow_dfield_getsampleresult(self._c, &P, &NP, &nq)
        cdef np.ndarray[np.double_t,ndim=2] x0 = np.zeros([points.shape[0], nd])
        cdef np.ndarray[np.double_t,ndim=2] P0 = np.zeros([points.shape[0], nq])
        cdef int i, j
        for i in range(Nx):
            for j in range(nd):
                x0[i,j] = x[3*i + j]
        for i in range(NP):
            for j in range(nq):
                P0[i,j] = P[nq*i + j]
        cow_dfield_setsamplecoords(self._c, NULL, 0, 3) # to save memory
        return x0, P0

    def index_global(self, ind):
        """
        Indexes the distributed array with the global index `ind`. Must be
        called on all ranks which own this array.
        """
        if len(ind) != self.domain.ndim:
            raise ValueError("wrong number of indices")
        i = [ind[n] if n < self.domain.ndim else 0 for n in range(3)]
        cdef double *P
        cdef int nq
        cow_dfield_sampleglobalind(self._c, i[0], i[1], i[2], &P, &nq)
        cdef np.ndarray[np.double_t,ndim=1] P0 = np.zeros(nq)
        cdef int m
        for m in range(nq):
            P0[m] = P[m]
        return P0

    cdef _apply_transform(self, args, cow_transform op, void *userdata=NULL):
        """
        Private method, removes redundant code in the two functions that follow.
        """
        cow_dfield_clearargs(self._c)
        for arg in args:
            cow_dfield_pusharg(self._c, (<DataField?>arg)._c)
        cow_dfield_settransform(self._c, op)
        cow_dfield_setuserdata(self._c, userdata)
        cow_dfield_transformexecute(self._c)
        return self

    def reduce_component(self, member):
        """
        Returns the min, max, and sum of the data member `member`, which may be
        int or str.
        """
        if type(member) is str:
            member = self.members.index(member)
        else:
            assert member < len(self.members)
        cdef np.ndarray[np.double_t,ndim=1] res = np.zeros(3)
        cow_dfield_clearargs(self._c)
        cow_dfield_settransform(self._c, cow_trans_component)
        cow_dfield_setuserdata(self._c, self._c)
        cow_dfield_setiparam(self._c, member)
        cow_dfield_reduce(self._c, <double*>res.data)
        return res
    
    def reduce_magnitude(self):
        """
        Returns the min, max, and sum of the data fields's vector magnitude.
        """
        cdef np.ndarray[np.double_t,ndim=1] res = np.zeros(3)
        cow_dfield_clearargs(self._c)
        cow_dfield_settransform(self._c, cow_trans_magnitude)
        cow_dfield_setuserdata(self._c, self._c)
        cow_dfield_reduce(self._c, <double*>res.data)
        return res

    def setflags_infnan(self):
        """
        Fills the data field's flags property (numpy integer array) with
        non-zero values wherever any of the data fields contain inf's or nan's.
        """
        cow_dfield_updateflaginfnan(self._c)

    def __getitem__(self, key):
        if type(key) is int:
            return self.value[..., key]
        else:
            return self.value[..., self.members.index(key)]

    def __setitem__(self, key, val):
        if type(key) is int:
            self.value[..., key] = val 
        else:
            self.value[..., self.members.index(key)] = val

    def __repr__(self):
        props = ["%s" % type(self),
                 "name: %s" % self.name,
                 "members: %s" % str(self.members),
                 "local shape: %s" % str(self.value[..., 0].shape),
                 "global shape: %s" % str(self.domain.global_shape),
                 "global start: %s" % str(self.domain.global_start),
                 "padding: %s" % self.domain.guard]
        return "{" + "\n\t".join(props) + "}"


cdef class ScalarField3d(DataField):
    def __init__(self, domain, members=("f",), name="scalarfield"):
        if len(members) != 1 or domain.ndim != 3:
            raise ValueError("bad argument list")
        super(ScalarField3d, self).__init__(domain, members, name)

    def power_spectrum(self, bins=128, spacing="linear", name=None):
        """
        Computes the spherically integrated power spectrum P(k) of the vector
        field, where P(k) = f(\vec{k}).
        """
        if name is None: name = self.name + "-pspec"
        cdef Histogram1d pspec = Histogram1d(0.0, 1.0, bins=bins,
                                             spacing=spacing, name=name,
                                             commit=False)
        cow_fft_pspecscafield(self._c, pspec._c)
        return pspec


cdef class VectorField3d(DataField):
    def __init__(self, domain, members=("fx","fy","fz"), name="vectorfield"):
        if len(members) != 3 or domain.ndim != 3:
            raise ValueError("bad argument list")
        super(VectorField3d, self).__init__(domain, members, name)

    def curl(self, name=None):
        """
        Takes the curl of the vector field using a 5-point stencil for the
        partial derivatives.
        """
        assert self.domain.guard >= 2
        if name is None: name = "del_cross_" + self.name
        cdef VectorField3d res = VectorField3d(self.domain, name=name)
        return res._apply_transform([self], cow_trans_rot5)

    def divergence(self, stencil="5point", name=None):
        """
        Takes the divergence of the vector field using a 5-point stencil for the
        partial derivatives if `stencil` is '5point', or the corner-valued 2nd
        order stencil of `stencil` is 'corner'.
        """
        if name is None: name = "del_dot_" + self.name
        cdef cow_transform op
        if stencil == "5point":
            assert self.domain.guard >= 2
            op = cow_trans_div5
        elif stencil == "corner":
            assert self.domain.guard >= 1
            op = cow_trans_divcorner
        else:
            raise ValueError("keyword 'stencil' must be one of ['5point', "
                             "'corner']")
        cdef ScalarField3d res = ScalarField3d(self.domain, name=name)
        return res._apply_transform([self], op)

    def solenoidal(self, name=None):
        """
        Returns the solenoidal (curl-like) part of the vector field's Helmholtz
        decomposition. Projection is done using the Fourier decomposition.
        """
        if name is None: name = self.name + "-solenoidal"
        cdef VectorField3d sol = VectorField3d(self.domain, name=name)
        sol.value[:] = self.value
        cow_fft_helmholtzdecomp(sol._c, COW_PROJECT_OUT_DIV)
        return sol

    def dilatational(self, name=None):
        """
        Returns the dilatational (div-like) part of the vector field's Helmholtz
        decomposition. Projection is done using the Fourier decomposition.
        """
        if name is None: name = self.name + "-dilatational"
        cdef VectorField3d dil = VectorField3d(self.domain, name=name)
        dil.value[:] = self.value
        cow_fft_helmholtzdecomp(dil._c, COW_PROJECT_OUT_CURL)
        return dil

    def power_spectrum(self, bins=128, spacing="linear", name=None):
        """
        Computes the spherically integrated power spectrum P(k) of the vector
        field, where P(k) = \vec{f}(\vec{k}) \cdot \vec{f}^*(\vec{k}).
        """
        if name is None: name = self.name + "-pspec"
        cdef Histogram1d pspec = Histogram1d(0.0, 1.0, bins=bins,
                                             spacing=spacing, name=name,
                                             commit=False)
        cow_fft_pspecvecfield(self._c, pspec._c)
        return pspec


cdef class Histogram1d(object):
    """
    Class that represents a 1 dimensional histogram.
    """
    _spacing = {"linear": COW_HIST_SPACING_LINEAR,
                "log": COW_HIST_SPACING_LOG}
    _binmode = {"counts": COW_HIST_BINMODE_COUNTS,
                "density": COW_HIST_BINMODE_DENSITY,
                "average": COW_HIST_BINMODE_AVERAGE}
    def __cinit__(self):
        self._c = cow_histogram_new()

    def __init__(self, x0, x1, bins=200, spacing="linear", binmode="counts",
                 name="histogram", domain=None, commit=True):
        cow_histogram_setlower(self._c, 0, x0)
        cow_histogram_setupper(self._c, 0, x1)
        cow_histogram_setnbins(self._c, 0, bins)
        cow_histogram_setbinmode(self._c, self._binmode[binmode])
        cow_histogram_setspacing(self._c, self._spacing[spacing])
        cow_histogram_setnickname(self._c, name)
        if domain:
            cow_histogram_setdomaincomm(self._c, (<DistributedDomain?>domain)._c)
        if commit:
            cow_histogram_commit(self._c)

    def __dealloc__(self):
        cow_histogram_del(self._c)

    property name:
        def __get__(self):
            return cow_histogram_getname(self._c)
        def __set__(self, value):
            cow_histogram_setnickname(self._c, value)

    @property
    def sealed(self):
        return bool(cow_histogram_getsealed(self._c))

    @property
    def counts(self):
        return cow_histogram_gettotalcounts(self._c)

    @property
    def binloc(self):
        """ Returns the bin centers """
        assert self.sealed
        cdef double *x
        cdef int N
        cow_histogram_getbinlocx(self._c, &x, &N)
        cdef np.ndarray[np.double_t,ndim=1] x0 = np.zeros(N)
        cdef int i
        for i in range(N):
            x0[i] = x[i]
        return x0

    @property
    def binval(self):
        """ Returns the present bin values """
        assert self.sealed
        cdef double *x
        cdef int N
        cow_histogram_getbinval1(self._c, &x, &N)
        cdef np.ndarray[np.double_t,ndim=1] x0 = np.zeros(N)
        cdef int i
        for i in range(N):
            x0[i] = x[i]
        return x0

    def add_sample(self, val, weight=1):
        """ Bins the data point with value `val` and weight `weight` """
        assert not self.sealed
        cow_histogram_addsample1(self._c, val, weight)

    def seal(self):
        """
        Locks out the addition of new samples. Also synchronizes across all
        participating ranks. Must be called before any function that gets
        histogram data or dumps it to file.
        """
        cow_histogram_seal(self._c)

    def dump(self, fname, gname="", format="guess"):
        """
        Writes the histogram to the HDF5 file `fname` under the group
        `gname`/self._name if `format` is 'hdf5', and an ascii table if
        'ascii'. By default `format` is guessed from the file extension.
        """
        assert self.sealed
        fname = str(fname)
        gname = str(gname)
        if format == "guess":
            if fname.endswith('.h5') or fname.endswith('.hdf5'):
                format = "hdf5" 
            elif fname.endswith('.dat') or fname.endswith('.txt'):
                format = "ascii"
            else:
                raise ValueError("file format could not be inferred from %s" %
                                 fname)
        if format == "hdf5":
            cow_histogram_dumphdf5(self._c, fname, gname)
        elif format == "ascii":
            cow_histogram_dumpascii(self._c, fname)
        else:
            raise ValueError("keyword 'format' must be one of ['hdf5', "
                             "'ascii']")


def dot_product(v, w):
    res = ScalarField3d(v.domain)
    res.name = v.name + "-dot-" + w.name
    return res._apply_transform([v, w], cow_trans_dot3)


def cross_product(v, w):
    res = VectorField3d(v.domain)
    res.name = v.name + "-cross-" + w.name
    return res._apply_transform([v, w], cow_trans_cross)


def fromfile(fname, group, guard=0, members=None, vec3d=False, downsample=0):
    """
    Looks in the HDF5 file `fname` for a serialized DataField named `group`, and
    creates a new data field from it if possible. All data sets in fname/group
    must have the same shape. If `members` is iterable, then only the data
    members it names will be read. If `vec3d` is True and there are exactly 3
    data members then a VectorField3d instance is returned.
    """
    h5f = h5py.File(fname, "r")
    grp = h5f[group]
    mem = [ ]
    shape = None
    for member in grp:
        if shape is None:
            shape = grp[member].shape # first time through
        if shape != grp[member].shape:
                raise ValueError("the group '%s' is not a DataField" % group)
        if not members or member in members:
            mem.append(member)
    if members:
        for m in members:
            if m not in mem:
                warnings.warn("data member '%s' was not found in file" % m)
    if downsample != 0:
        sprime = [int(np.ceil(float(s) / downsample)) for s in shape]
        domain = DistributedDomain(sprime, guard=guard)
        if domain.cart_size != 1:
            raise ValueError("down-sampled loading only supported for serial"
                             " jobs right now")
    else:
        domain = DistributedDomain(shape, guard=guard)
    if vec3d:
        assert len(mem) == 3
        dfield = VectorField3d(domain, mem, name=group)
    else:
        dfield = DataField(domain, mem, name=group)
    if downsample == 0:
        h5f.close() # make sure to close before opening another instance
        dfield.read(fname)
    else:
        s = downsample
        for n, m in enumerate(mem):
            if domain.ndim == 1:
                dfield.interior[...,n] = h5f[group][m][::s]
            elif domain.ndim == 2:
                dfield.interior[...,n] = h5f[group][m][::s,::s]
            elif domain.ndim == 3:
                dfield.interior[...,n] = h5f[group][m][::s,::s,::s]
        h5f.close()
        dfield.sync_guard()
    return dfield
