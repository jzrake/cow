
import numpy as np
import cowpy

def testsamp():
    domain = cowpy.DistributedDomain([10,10,10], guard=3)
    dfield = cowpy.DataField(domain, ["vx", "vy", "vz"])
    np.random.seed(domain.rank)
    sampx = np.random.rand(10,3)
    dfield["vx"] = 2.1
    dfield["vy"] = 2.2
    dfield["vz"] = 2.3
    x, P = dfield.sample(sampx)
    assert abs(P[:,0] - 2.1).all() < 1e-16
    assert abs(P[:,1] - 2.2).all() < 1e-16
    assert abs(P[:,2] - 2.3).all() < 1e-16
    assert len(x) == len(sampx)

def testio():
    domain = cowpy.DistributedDomain([12,12,16], guard=3)
    dfield = cowpy.DataField(domain, ["vx", "vy", "vz"], name="vel")
    dfield["vx"] = 2.1
    dfield["vy"] = 2.2
    dfield["vz"] = 2.3
    dfield.dump("iotest.h5")
    dfield.value[:] = 0.0
    dfield.read("iotest.h5")
    assert abs(dfield["vx"] - 2.1).all() < 1e-16
    assert abs(dfield["vy"] - 2.2).all() < 1e-16
    assert abs(dfield["vz"] - 2.3).all() < 1e-16

def testhist():
    hist = cowpy.Histogram1d(0, 1, N=10, binmode="counts")
    for n in range(1000):
        hist.add_sample(np.random.rand())
    print hist.binloc
    assert abs(hist.binval.sum() - 1000) < 1e-16

def testcb():
    domain = cowpy.DistributedDomain([10,10,10], guard=2)
    B = cowpy.VectorField3d(domain, name="B")
    B[0] = 1.0
    B[1] = 1.1
    B[2] = 1.2
    J = B.curl()
    M1 = B.divergence(method="corner")
    M2 = B.divergence(method="5point")
    print J
    print M1

if __name__ == "__main__":
    testhist()
    testsamp()
    testio()
    testcb()
