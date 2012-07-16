
import os
import atexit
import numpy as np
from _cow import *


def cow_exit():
    print "finalizing cow..."
    cow_finalize()

KILOBYTES = 1 << 10
MEGABYTES = 1 << 20
modes = 0
hdf5_collective = os.getenv("COW_HDF5_COLLECTIVE", 0)
modes |= (COW_NOREOPEN_STDOUT if os.getenv("COW_NOREOPEN_STDOUT", 0) else 0)
modes |= (COW_DISABLE_MPI if os.getenv("COW_DISABLE_MPI", 0) else 0)

atexit.register(cow_exit)
cow_init(0, None, modes)


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
    def rank(self):
        return cow_domain_getcartrank(self._cdomain)
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
        assert type(name) is str
        self._cdfield = cow_dfield_new(domain._cdomain, name)
        for m in members:
            assert(type(m) is str)
            cow_dfield_addmember(self._cdfield, m)
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
    def value(self):
        return self._buf
    @property
    def domain(self):
        return self._domain

    def dump(self, fname):
        assert type(fname) is str
        cow_dfield_write(self._cdfield, fname)

    def read(self, fname):
        assert type(fname) is str
        cow_dfield_read(self._cdfield, fname)

    def sample(self, points):
        cow_dfield_setsamplecoords(self._cdfield, points)
        cow_dfield_setsamplemode(self._cdfield, COW_SAMPLE_LINEAR)
        cow_dfield_sampleexecute(self._cdfield)
        x = cow_dfield_getsamplecoords(self._cdfield).copy()
        P = cow_dfield_getsampleresult(self._cdfield).copy()
        return  x, P

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
    def __init__(self, domain, members=["f"], name="scalarfield"):
        assert len(members) == 1
        assert domain._nd == 3
        super(ScalarField3d, self).__init__(domain, members, name)


class VectorField3d(DataField):
    def __init__(self, domain, members=["fx","fy","fz"], name="vectorfield"):
        assert len(members) == 3
        assert domain._nd == 3
        super(VectorField3d, self).__init__(domain, members, name)

    def curl(self):
        assert self.domain.guard >= 2
        res = VectorField3d(self._domain, name="del_cross_"+self._name)
        cow_dfield_transf1(res._cdfield, self._cdfield, cow_trans_rot5, None)
        return res

    def divergence(self, method="5point"):
        trans = cow_dfield_transf1
        if method == "5point":
            assert self.domain.guard >= 2
            res = ScalarField3d(self._domain, name="del_dot_"+self._name)
            trans(res._cdfield, self._cdfield, cow_trans_div5, None)
            return res
        elif method == "corner":
            assert self.domain.guard >= 1
            res = ScalarField3d(self._domain, name="del_dot_"+self._name)
            trans(res._cdfield, self._cdfield, cow_trans_divcorner, None)
            return res
        else:
            raise ValueError("keyword 'method' must be one of ['5point', "
                             "'corner']")

class Histogram1d(object):
    """
    Class that represents a 1 dimensional histogram.
    """
    _spacing = {"linear": COW_HIST_SPACING_LINEAR,
                "log": COW_HIST_SPACING_LOG}
    _binmode = {"counts": COW_HIST_BINMODE_COUNTS,
                "density": COW_HIST_BINMODE_DENSITY,
                "average": COW_HIST_BINMODE_AVERAGE}

    def __init__(self, x0, x1, N=200, spacing="linear", binmode="counts",
                 name="histogram"):
        self._chist = cow_histogram_new()
        self._name = name
        cow_histogram_setlower(self._chist, 0, x0)
        cow_histogram_setupper(self._chist, 0, x1)
        cow_histogram_setnbins(self._chist, 0, N)
        cow_histogram_setbinmode(self._chist, self._binmode[binmode])
        cow_histogram_setspacing(self._chist, self._spacing[spacing])
        cow_histogram_setnickname(self._chist, name)
        cow_histogram_commit(self._chist)

    def __del__(self):
        cow_histogram_del(self._chist)

    def add_sample(self, val, weight=1):
        """ Bins the data point with value `val` and weight `weight` """
        cow_histogram_addsample1(self._chist, val, weight)

    def dump(self, fname, gname="", format=None):
        """
        Writes the histogram to the HDF5 file `fname` under the group
        `gname`/self._name if `format` is 'hdf5', and an ascii table if
        'ascii'. By default `format` to guessed from the file extension.
        """
        assert type(fname) is str
        assert type(gname) is str
        if not format:
            if fname.endswith('.h5'):
                format = "hdf5" 
            elif fname.endswith('.dat') or fname.endswith('.txt'):
                format = "ascii" 
        if format == "hdf5":
            cow_histogram_dumphdf5(self._chist, fname, gname)
        elif format == "ascii":
            cow_histogram_dumpascii(self._chist, fname)
        else:
            raise ValueError("keyword 'format' must be one of ['hdf5', "
                             "'ascii']")

    @property
    def binloc(self):
        """ Returns the bin centers """
        return cow_histogram_getbinlocx(self._chist).copy()

    @property
    def binval(self):
        """ Returns the present bin values """
        return cow_histogram_getbinval1(self._chist).copy()
