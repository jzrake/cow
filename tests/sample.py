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
    """
    Demonstrates that the samples coorindates that go out are the same ones that
    come back, as long as they are re-sorted.
    """
    ng = 3
    nx = 32
    domain = cowpy.DistributedDomain([nx,nx,nx], guard=ng)
    dfield = cowpy.VectorField3d(domain)
    ishape = dfield.interior.shape
    oshape = dfield.value.shape
    np.random.seed(domain.cart_rank)

    for i,j,k in itertools.product(*(range(n) for n in oshape[0:3])):
        dfield.value[i,j,k] = domain.coordinate([i,j,k])
    dfield.sync_guard()

    # Syncing guard zones will throw off the sampling if they're near the
    # boundary, so choose them suffieciently far away. This is not a problem
    # when the data is periodic.
    nsamp = 100
    sampx = 0.125 + 0.25*np.random.rand(nsamp, 3)

    # Indexing trick to sort array by first index:
    # http://stackoverflow.com/questions/2828059/sorting-arrays-in-numpy-by-column
    x, P = dfield.sample_global(sampx)
    xi = sampx[sampx[:,0].argsort(),:] # x in
    xo = x[x[:,0].argsort(),:] # x out
    sP = P[P[:,0].argsort(),:]

    # The samples that go out are the same that come back
    assert (abs(xo - xi) < 1e-14).all()

    # The samples are correctly matched with their coordinates
    assert (abs(x - P) < 1e-14).all()


if __name__ == "__main__":
    test2()
