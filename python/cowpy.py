
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
from _cow import *


KILOBYTES = 1 << 10
MEGABYTES = 1 << 20
modes = 0
hdf5_collective = os.getenv("COW_HDF5_COLLECTIVE", 0)
modes |= (COW_NOREOPEN_STDOUT if int(os.getenv("COW_NOREOPEN_STDOUT", 0)) else 0)
modes |= (COW_DISABLE_MPI if int(os.getenv("COW_DISABLE_MPI", 0)) else 0)
modes |= (COW_DISABLE_MPI if '-s' in sys.argv else 0)

def cow_exit():
    print "finalizing cow..."
    cow_finalize()

atexit.register(cow_exit)
cow_init(0, None, modes)

mpirunning = cow_mpirunning

class DistributedDomain(object):
    def __init__(self, G_ntot, guard=0, x0=None, x1=None):
        self._nd = len(G_ntot)
        assert self._nd <= 3
        self._cdomain = cow_domain_new()
        self._G_ntot = tuple(G_ntot)
        for n, ni in enumerate(G_ntot):
            cow_domain_setsize(self._cdomain, n, ni)
        cow_domain_setndim(self._cdomain, self._nd)
        cow_domain_setguard(self._cdomain, guard)
        cow_domain_commit(self._cdomain)
        # Set up the IO scheme
        cow_domain_setchunk(self._cdomain, 1)
        cow_domain_setcollective(self._cdomain, hdf5_collective)
        cow_domain_setalign(self._cdomain, 4*KILOBYTES, 4*MEGABYTES)

    def __del__(self):
        cow_domain_del(self._cdomain)

    @property
    def cart_rank(self):
        return cow_domain_getcartrank(self._cdomain)

    @property
    def cart_size(self):
        return cow_domain_getcartsize(self._cdomain)

    @property
    def global_start(self):
        return tuple(cow_domain_getglobalstartindex(self._cdomain, n)
                     for n in range(self._nd))
    @property
    def guard(self):
        return cow_domain_getguard(self._cdomain)

    def coordinate(self, ind):
        assert len(ind) == self._nd
        return [cow_domain_positionatindex(self._cdomain, n, i)
                for n, i in enumerate(ind)]


class DataField(object):
    def __init__(self, domain, members, name="datafield"):
        self._cdfield = cow_dfield_new(domain._cdomain, name)
        for m in members:
            cow_dfield_addmember(self._cdfield, str(m))
        nd = cow_domain_getndim(domain._cdomain)
        dims = [ ]
        for i in range(nd):
            dims.append(cow_domain_getnumlocalzonesincguard(domain._cdomain, i))
        dims.append(len(members))
        self._buf = np.zeros(dims)
        if nd == 1:
            setarray1(self._cdfield, self._buf)
        if nd == 2:
            setarray2(self._cdfield, self._buf)
        if nd == 3:
            setarray3(self._cdfield, self._buf)
        cow_dfield_commit(self._cdfield)
        self._domain = domain
        self._members = tuple(members)
        self._lookup = dict([(m, n) for n,m in enumerate(members)])
        self._name = name

    def __del__(self):
        cow_dfield_del(self._cdfield)

    @property
    def name(self):
        return cow_dfield_getname(self._cdfield)

    @property
    def value(self):
        """
        Returns the numpy.ndarray object serving as the buffer for the
        underlying C API. The last dimension contains the data members, so the
        result is an (N+1)-d array for an N-d domain. The first and last zones
        in the array are padding, or 'guard' zones.
        """
        return self._buf

    @property
    def interior(self):
        """
        Returns a view of the buffer object which refers only to the interior
        portion, excluding the guard zones.
        """
        ng = self.domain.guard
        if ng == 0:
            return self._buf
        elif self.domain._nd == 1:
            return self._buf[ng:-ng, :]
        elif self.domain._nd == 2:
            return self._buf[ng:-ng, ng:-ng, :]
        elif self.domain._nd == 3:
            return self._buf[ng:-ng, ng:-ng, ng:-ng, :]

    @property
    def domain(self):
        return self._domain

    def dump(self, fname):
        cow_dfield_write(self._cdfield, str(fname))

    def read(self, fname):
        cow_dfield_read(self._cdfield, str(fname))

    def sample_global(self, points):
        """
        Samples the distributed domain at the provided list of physical
        coordinates, `points` which has the shape (N x nd). Returns those
        coordinates, but permuted arbitrarily, and their associated field
        values. Values are multi-linearly interpolated from cell centers. This
        function must be called by all ranks which own this array, although the
        size of `points` may vary between ranks, and is allowed to be zero.
        """
        assert points.shape[1] == self.domain._nd
        p3d = np.zeros([points.shape[0], 3])
        for n in range(self.domain._nd):
            p3d[:,n] = points[:,n]
        cow_dfield_setsamplecoords(self._cdfield, p3d)
        cow_dfield_setsamplemode(self._cdfield, COW_SAMPLE_LINEAR)
        cow_dfield_sampleexecute(self._cdfield)
        x = cow_dfield_getsamplecoords(self._cdfield).copy()
        P = cow_dfield_getsampleresult(self._cdfield).copy()
        # To save memory, reset the sample coordinate buffer
        cow_dfield_setsamplecoords(self._cdfield, np.zeros([0,3]))
        return  x, P

    def index_global(self, ind):
        """
        Indexes the distributed array with the global index `ind`. Must be
        called on all ranks which own this array.
        """
        assert len(ind) == self.domain._nd
        i = [ind[n] if n < self.domain._nd else 0 for n in range(3)]
        return cow_dfield_sampleglobalind(self._cdfield, i[0], i[1], i[2])

    def _apply_transform(self, args, op, userdata=None):
        cow_dfield_clearargs(self._cdfield)
        for arg in args:
            cow_dfield_pusharg(self._cdfield, arg._cdfield)
        cow_dfield_settransform(self._cdfield, op)
        cow_dfield_setuserdata(self._cdfield, userdata)
        cow_dfield_transformexecute(self._cdfield)
        return self

    def reduce_component(self, member):
        """
        Returns the min, max, and sum of the data member `member`, which may be
        int or str.
        """
        if type(member) is str:
            member = self._lookup[member]
        else:
            assert member < len(self._members)
        cow_dfield_clearargs(self._cdfield)
        cow_dfield_settransform(self._cdfield, cow_trans_component)
        cow_dfield_setuserdata(self._cdfield, self._cdfield)
        cow_dfield_setiparam(self._cdfield, member)
        return cow_dfield_reduce(self._cdfield)

    def reduce_magnitude(self):
        """
        Returns the min, max, and sum of the data fields's vector magnitude.
        """
        cow_dfield_clearargs(self._cdfield)
        cow_dfield_settransform(self._cdfield, cow_trans_magnitude)
        cow_dfield_setuserdata(self._cdfield, self._cdfield)
        return cow_dfield_reduce(self._cdfield)

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
                 "members: %s" % str(self._members),
                 "local shape: %s" % str(self.value[..., 0].shape),
                 "global shape: %s" % str(self._domain._G_ntot),
                 "global start: %s" % str(self._domain.global_start),
                 "padding: %s" % self.domain.guard]
        return "{" + "\n\t".join(props) + "}"


class ScalarField3d(DataField):
    def __init__(self, domain, members=("f",), name="scalarfield"):
        assert len(members) == 1
        assert domain._nd == 3
        super(ScalarField3d, self).__init__(domain, members, name)


class VectorField3d(DataField):
    def __init__(self, domain, members=("fx","fy","fz"), name="vectorfield"):
        assert len(members) == 3
        assert domain._nd == 3
        super(VectorField3d, self).__init__(domain, members, name)

    def curl(self, name=None):
        """
        Takes the curl of the vector field using a 5-point stencil for the
        partial derivatives.
        """
        assert self.domain.guard >= 2
        if name is None: name = "del_cross_" + self._name
        res = VectorField3d(self._domain, name=name)
        return res._apply_transform([self], cow_trans_rot5)

    def divergence(self, stencil="5point", name=None):
        """
        Takes the divergence of the vector field using a 5-point stencil for the
        partial derivatives if `stencil` is '5point', or the corner-valued 2nd
        order stencil of `stencil` is 'corner'.
        """
        if name is None: name = "del_dot_" + self._name
        if stencil == "5point":
            assert self.domain.guard >= 2
            op = cow_trans_div5
        elif stencil == "corner":
            assert self.domain.guard >= 1
            op = cow_trans_divcorner
        else:
            raise ValueError("keyword 'stencil' must be one of ['5point', "
                             "'corner']")
        res = ScalarField3d(self._domain, name=name)
        return res._apply_transform([self], op)

    def solenoidal(self, name=None):
        """
        Returns the solenoidal (curl-like) part of the vector field's Helmholtz
        decomposition. Projection is done using the Fourier decomposition.
        """
        if name is None: name = self._name + "-solenoidal"
        sol = VectorField3d(self._domain, members=self._members, name=name)
        sol.value[:] = self.value
        cow_fft_helmholtzdecomp(sol._cdfield, COW_PROJECT_OUT_DIV)
        return sol

    def dilatational(self, name=None):
        """
        Returns the dilatational (div-like) part of the vector field's Helmholtz
        decomposition. Projection is done using the Fourier decomposition.
        """
        if name is None: name = self._name + "-dilatational"
        div = VectorField3d(self._domain, members=self._members, name=name)
        div.value[:] = self.value
        cow_fft_helmholtzdecomp(div._cdfield, COW_PROJECT_OUT_CURL)
        return div

    def power_spectrum(self, bins=200, spacing="linear", name=None):
        """
        Computes the spherically integrated power spectrum P(k) of the vector
        field, where P(k) = \vec{f}(\vec{k}) \cdot \vec{f}^*(\vec{k}).
        """
        if name is None: name = self._name + "-pspec"
        pspec = Histogram1d(0.0, 1.0, bins=bins, spacing=spacing,
                            name=name, commit=False)
        cow_fft_pspecvecfield(self._cdfield, pspec._chist)
        return pspec


class Histogram1d(object):
    """
    Class that represents a 1 dimensional histogram.
    """
    _spacing = {"linear": COW_HIST_SPACING_LINEAR,
                "log": COW_HIST_SPACING_LOG}
    _binmode = {"counts": COW_HIST_BINMODE_COUNTS,
                "density": COW_HIST_BINMODE_DENSITY,
                "average": COW_HIST_BINMODE_AVERAGE}

    def __init__(self, x0, x1, bins=200, spacing="linear", binmode="counts",
                 name="histogram", domain=None, commit=True):
        self._chist = cow_histogram_new()
        self._name = name
        cow_histogram_setlower(self._chist, 0, x0)
        cow_histogram_setupper(self._chist, 0, x1)
        cow_histogram_setnbins(self._chist, 0, bins)
        cow_histogram_setbinmode(self._chist, self._binmode[binmode])
        cow_histogram_setspacing(self._chist, self._spacing[spacing])
        cow_histogram_setnickname(self._chist, name)
        if domain:
            cow_histogram_setdomaincomm(self._chist, domain._cdomain)
        if commit:
            cow_histogram_commit(self._chist)

    def __del__(self):
        cow_histogram_del(self._chist)

    @property
    def sealed(self):
        return bool(cow_histogram_getsealed(self._chist))

    @property
    def binloc(self):
        """ Returns the bin centers """
        assert self.sealed
        return cow_histogram_getbinlocx(self._chist).copy()

    @property
    def binval(self):
        """ Returns the present bin values """
        assert self.sealed
        return cow_histogram_getbinval1(self._chist).copy()

    def add_sample(self, val, weight=1):
        """ Bins the data point with value `val` and weight `weight` """
        assert not self.sealed
        cow_histogram_addsample1(self._chist, val, weight)

    def seal(self):
        """
        Locks out the addition of new samples. Also synchronizes across all
        participating ranks. Must be called before any function that gets
        histogram data or dumps it to file.
        """
        cow_histogram_seal(self._chist)

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
            cow_histogram_dumphdf5(self._chist, fname, gname)
        elif format == "ascii":
            cow_histogram_dumpascii(self._chist, fname)
        else:
            raise ValueError("keyword 'format' must be one of ['hdf5', "
                             "'ascii']")

    def __setattr__(self, key, value):
        if key == "name":
            cow_histogram_setnickname(self._chist, value)
        else:
            object.__setattr__(self, key, value)


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
    for member in grp:
        try:
            if shape != grp[member].shape:
                raise AttributeError # happens if grp[member] is not a dataset
        except AttributeError:
                raise ValueError("the group '%s' is not a DataField" % group)
        except NameError:
            shape = grp[member].shape # first time through
        if not members or member in members:
            mem.append(member)
    if members:
        for m in members:
            if m not in mem:
                warnings.warn("data member '%s' was not found in file" % m)
    if downsample != 0:
        sprime = [s / downsample for s in shape]
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
    if downsample:
        s = downsample
        for m in mem:
            if domain._nd == 1:
                dfield[m][...] = h5f[group][m][::s]
            elif domain._nd == 2:
                dfield[m][...] = h5f[group][m][::s,::s]
            elif domain._nd == 3:
                dfield[m][...] = h5f[group][m][::s,::s,::s]
        del s
    else:
        dfield.read(fname)

    h5f.close()
    return dfield
