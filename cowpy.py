
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


class UnigridDatafield(object):
    def __init__(self, domain, members, name="dfield"):
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

    def __del__(self):
        print "deleting dfield"
        cow_dfield_del(self._cdfield)

    @property
    def value(self):
        return self._buf

    def dump(self, fname):
        assert type(fname) is str
        cow_dfield_write(self._cdfield, fname)

    def sample(self, points):
        cow_dfield_setsamplecoords(self._cdfield, points)
        cow_dfield_setsamplemode(self._cdfield, COW_SAMPLE_LINEAR)
        cow_dfield_sampleexecute(self._cdfield)
        x = cow_dfield_getsamplecoords(self._cdfield).copy()
        P = cow_dfield_getsampleresult(self._cdfield).copy()
        return  x, P

class UnigridDomain(object):
    def __init__(self, G_ntot, guard=0, x0=None, x1=None):
        self._nd = len(G_ntot)
        assert self._nd <= 3
        self._cdomain = cow_domain_new()
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
        print "deleting domain"
        cow_domain_del(self._cdomain)

    @property
    def rank(self):
        return cow_domain_getcartrank(self._cdomain)

    @property
    def global_start(self):
        return [cow_domain_getglobalstartindex(self._cdomain, n)
                for n in range(self._nd)]

    def coordinate(self, *args):
        assert len(args) == self._nd
        return [cow_domain_positionatindex(self._cdomain, n, i)
                for n, i in enumerate(args)]
