
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


class UnigridDatafield:
    def __init__(self, domain, members, name="dfield"):
        assert(type(name) is str)
        self._cdfield = cow_dfield_new(domain._cdomain, name)
        for m in members:
            assert(type(m) is str)
            cow_dfield_addmember(self._cdfield, m)
        cow_dfield_commit(self._cdfield)
        self._domain = domain
        self.value = cow_dfield_getarray(self._cdfield)

    def __del__(self):
        print "deleting dfield"
        cow_dfield_del(self._cdfield)


class UnigridDomain:
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
    print dfield.value.max()


if __name__ == "__main__":
    test()
