
import os
import numpy as np
import cowpy

N = int(os.getenv("COW_TESTPARALLELIO_RES", 64))
fname = os.getenv("COW_TESTPARALLELIO_OUTFILE", "data/testpario.h5")

def test1():
    domain = cowpy.DistributedDomain([N, N, N], guard=3)
    dfield = cowpy.DataField(domain, ["vx"], name="vel1")
    dfield.dump(fname)

def test2():
    domain = cowpy.DistributedDomain([N, N, N], guard=3)
    dfield = cowpy.DataField(domain, ["vx", "vy"], name="vel2")
    dfield.dump(fname)

def test3():
    domain = cowpy.DistributedDomain([N, N, N], guard=3)
    dfield = cowpy.DataField(domain, ["vx", "vy", "vz"], name="vel3")
    dfield.dump(fname)

def test4():
    domain = cowpy.DistributedDomain([N, N, N], guard=3)
    dfield = cowpy.VectorField3d(domain, name="vel4")
    dfield.dump(fname)

def test5():
    domain = cowpy.DistributedDomain([N, N, N], guard=3)
    dfield = cowpy.VectorField3d(domain, name="vel4")
    div = dfield.divergence(name="div")
    rot = dfield.curl(name="rot")
    dfield.dump(fname)
    div.dump(fname)
    rot.dump(fname)


if __name__ == "__main__":
    test1()
    test2()
    test3()
    test4()
    test5()
