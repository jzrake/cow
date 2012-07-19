#!/usr/bin/env python

import cowpy
import numpy as np
import itertools

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
    domain = cowpy.DistributedDomain([16,16,16], guard=ng)
    dfield = cowpy.VectorField3d(domain)
    ishape = dfield.interior.shape
    oshape = dfield.value.shape

    for i,j,k in itertools.product(*(range(n) for n in oshape[0:3])):
        dfield.value[i,j,k] = domain.coordinate([i,j,k])
    dfield.sync_guard()

    sampx = [domain.coordinate([i+ng, j+ng, k+ng]) for i,j,k in
             itertools.product(*(range(n) for n in ishape[0:3]))]
    x, P = dfield.sample_global(sampx)
    print x[0], P[0]
    print x[1], P[1]
    print x[2], P[2]


if __name__ == "__main__":
    test2()
