#!/usr/bin/env python

import sys
import itertools
import numpy as np
import cowpy

# General convention is that local indices include guard zones, and global ones
# do not.

def test1():
    ng = 2
    domain = cowpy.DistributedDomain([16,16,16], guard=ng)
    dfield = cowpy.VectorField3d(domain)
    ishape = dfield.interior.shape
    oshape = dfield.value.shape

    i0, j0, k0 = domain.global_start

    # domain.coordinate returns the position at the requested index- which is
    # relative to the local subgrid, and has guard zones added to it, i.e. i=ng
    # means x=dx/2.
    for i,j,k in itertools.product(*(range(n) for n in oshape[0:3])):
        dfield.value[i,j,k] = domain.coordinate([i,j,k])
    dfield.sync_guard()

    errors = 0
    for i,j,k in itertools.product(*(range(n) for n in ishape[0:3])):
        x = np.array(domain.coordinate([i+ng, j+ng, k+ng]))
        y = dfield.index_global([i+i0, j+j0, k+k0])
        try:
            assert (x == y).all()
        except AssertionError:
            print domain.cart_rank, x, y
            errors += 1
    assert errors == 0


def test2():
    ng = 1
    nx = 32
    domain = cowpy.DistributedDomain([nx,nx,nx], guard=ng)
    dfield = cowpy.VectorField3d(domain)
    ishape = dfield.interior.shape
    oshape = dfield.value.shape
    np.random.seed(domain.cart_rank)

    for i,j,k in itertools.product(*(range(n) for n in oshape[0:3])):
        dfield.value[i,j,k] = domain.coordinate([i,j,k])

    domain.sequential(lambda: sys.stdout.write("before: %s\n" % dfield.value[-1,-1,-1]))
    dfield.sync_guard()
    domain.sequential(lambda: sys.stdout.write("after:  %s\n" % dfield.value[-1,-1,-1]))

    nsamp = 10000
    npair = 100000
    sampx = np.random.rand(nsamp, 3)

    x, P = dfield.sample_global(sampx)
    histx = cowpy.Histogram1d(0.0, 1.5, bins=36, binmode="counts")
    histP = cowpy.Histogram1d(0.0, 1.5, bins=36, binmode="counts")

    for n in range(npair):
        i0 = np.random.randint(0, len(P))
        i1 = np.random.randint(0, len(P))
        histx.add_sample(np.dot(x[i1] - x[i0], x[i1] - x[i0])**0.5)
        histP.add_sample(np.dot(P[i1] - P[i0], P[i1] - P[i0])**0.5)
    histx.seal()
    histP.seal()

    if domain.cart_rank == 0:
        import matplotlib.pyplot as plt
        plt.plot(histx.binloc, histx.binval, label="x")
        plt.plot(histP.binloc, histP.binval, label="P")
        plt.show()


if __name__ == "__main__":
    test1()
