
from _cow import *
import os, sys
import atexit
import numpy as np
from sys import getrefcount


def cow_exit():
    print "finalizing cow..."
    cow_finalize()

atexit.register(cow_exit)
cow_init()


class UnigridDatafield(object):
    def __init__(self, domain, members, name="dfield"):
        assert(type(name) is str)
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
        assert(type(fname) is str)
        cow_dfield_write(self._cdfield, fname)


class UnigridDomain(object):
    def __init__(self, G_ntot, guard=0, x0=None, x1=None):
        self._nd = len(G_ntot)
        assert(self._nd <= 3)
        self._cdomain = cow_domain_new()
        for n, ni in enumerate(G_ntot):
            cow_domain_setsize(self._cdomain, n, ni)
        cow_domain_setndim(self._cdomain, self._nd)
        cow_domain_setguard(self._cdomain, guard)
        cow_domain_commit(self._cdomain)

    def __del__(self):
        print "deleting domain"
        cow_domain_del(self._cdomain)


def test():
    domain = UnigridDomain([10,10,10], guard=3)
    dfield = UnigridDatafield(domain, ["vx", "vy", "vz"])
    dfield.value[4:5,:,:,0] = 2.2
    dfield.dump("pyout.h5")
    print dfield.value.max()


if __name__ == "__main__":
    test()


