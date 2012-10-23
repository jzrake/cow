
import sys
import numpy as np
import cowpy

"""
All tests here are written assuming one MPI process.
"""

def testfromfile():
    field = cowpy.fromfile("test.h5", group="vel")
    print field
    sparsefield = cowpy.fromfile("test.h5", group="vel", downsample=2)
    assert sparsefield.value.shape == (12/2, 12/2, 16/2, 3)
    print sparsefield
    # Test with a strange downsampling factor and with guard cells
    sparsefield = cowpy.fromfile("test.h5", group="vel", downsample=3, guard=3)
    print sparsefield
    assert sparsefield.value.shape == (4+2*3, 4+2*3, 6+2*3, 3)

def testsamp():
    domain = cowpy.DistributedDomain([10,10,10], guard=3)
    dfield = cowpy.DataField(domain, ["vx", "vy", "vz"])
    np.random.seed(domain.cart_rank)
    sampx = np.random.rand(10,3)
    dfield["vx"] = 2.1
    dfield["vy"] = 2.2
    dfield["vz"] = 2.3
    x, P = dfield.sample_global(sampx)
    assert abs(P[:,0] - 2.1).all() < 1e-16
    assert abs(P[:,1] - 2.2).all() < 1e-16
    assert abs(P[:,2] - 2.3).all() < 1e-16
    assert len(x) == len(sampx)

def testglobind():
    domain = cowpy.DistributedDomain([10,10,10], guard=3)
    dfield = cowpy.DataField(domain, ["vx", "vy", "vz"])
    dfield.interior[:,:,:,:] = np.mgrid[0:10,0:10,0:10].transpose([1,2,3,0])
    x, y, z = dfield.index_global([2, 3, 4])
    assert abs(x - 2) < 1e-16
    assert abs(y - 3) < 1e-16
    assert abs(z - 4) < 1e-16

def testio():
    domain = cowpy.DistributedDomain([12,12,16], guard=3)
    dfield = cowpy.DataField(domain, ["vx", "vy", "vz"], name="vel")
    dfield["vx"] = 2.1
    dfield["vy"] = 2.2
    dfield["vz"] = 2.3
    dfield.dump("test.h5")
    dfield.value[:] = 0.0
    dfield.read("test.h5")
    assert abs(dfield["vx"] - 2.1).all() < 1e-16
    assert abs(dfield["vy"] - 2.2).all() < 1e-16
    assert abs(dfield["vz"] - 2.3).all() < 1e-16

def testhist():
    domain = cowpy.DistributedDomain([12,12,16], guard=3)
    hist = cowpy.Histogram1d(0, 1, bins=10, binmode="counts", name="myhist",
                             domain=domain)
    for n in range(1000):
        hist.add_sample(np.random.rand())
    hist.seal()
    assert abs(hist.binval.sum() - 1000) < 1e-16
    hist.name = "myrealhist"
    hist.dump("hist.dat")
    hist.dump("test.h5", group="G1/G2")

def testmaghist():
    domain = cowpy.DistributedDomain([12,12,16], guard=3)
    dfield = cowpy.VectorField3d(domain)
    dfield.value[...] = np.random.rand(*dfield.value.shape)
    min, max, sum = dfield.reduce_magnitude()
    maghist = cowpy.Histogram1d(min, max, bins=100, binmode="counts", name="myhist",
                                domain=domain)    
    vx, vy, vz = dfield[0], dfield[1], dfield[2]
    mag = (vx**2 + vy**2 + vz**2)**0.5
    for a in mag.flatten(): maghist.add_sample(a)
    maghist.seal()
    print maghist.binval

def testhelm():
    domain = cowpy.DistributedDomain([10,10,10], guard=2)
    A = cowpy.VectorField3d(domain, name="B")
    Asol = A.solenoidal()
    Adil = A.dilatational()
    pspec = Adil.power_spectrum(bins=100, spacing="log")
    Asol.dump("test.h5")
    Adil.dump("test.h5")
    pspec.dump("test.h5")

def testcb():
    domain = cowpy.DistributedDomain([10,10,10], guard=2)
    B = cowpy.VectorField3d(domain, name="B")
    B[0] = 1.0
    B[1] = 1.1
    B[2] = 1.2
    J = B.curl()
    M1 = B.divergence(stencil="corner")
    M2 = B.divergence(stencil="5point")
    assert J.name == "del_cross_B"
    assert M1.name == "del_dot_B"

def testreduce():
    domain = cowpy.DistributedDomain([10,10,10], guard=2)
    J = cowpy.VectorField3d(domain, name="J")
    J[0] = 0.0
    J[0][2,0,2] = 20.0 # guard cells should not be included in the reduction
    J[0][2,12,2] = -20.0
    J[0][2,2,2] = 10.0
    J[0][2,2,3] =-10.0
    min, max, sum = J.reduce_component("fx")
    assert abs(min + 10.0) < 1e-16
    assert abs(max - 10.0) < 1e-16
    assert abs(sum - 0.0) < 1e-16
    min, max, sum = J.reduce_magnitude()
    assert abs(max - 10.0) < 1e-16
    assert abs(min - 0.0) < 1e-16

def testbadsetup():
    try:
        e = None
        domain = cowpy.DistributedDomain([6,6,6,5], guard=2)
    except ValueError as e:
        print e
    assert isinstance(e, ValueError)

    try:
        e = None
        domain = cowpy.DistributedDomain([6,6,6], guard=2)
        vfield = cowpy.VectorField3d(domain, name=3)
    except TypeError as e:
        print e
    assert isinstance(e, TypeError)

    try:
        e = None
        domain = cowpy.DistributedDomain([6,6,6], guard=2)
        vfield = cowpy.VectorField3d(domain, members=["gx", "gy"])
    except ValueError as e:
        print e
    assert isinstance(e, ValueError)

def testuserbuffer():
    buffer = np.zeros([10,10,10,1], order='C')
    domain = cowpy.DistributedDomain([6,6,6], guard=2)
    dfield = cowpy.DataField(domain, members=['f'], buffer=buffer)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        print "running test", sys.argv[1]
        globals()["test" + sys.argv[1]]()
    else:
        print "running all tests..."
        testhist()
        testmaghist()
        testsamp()
        testglobind()
        testio()
        testcb()
        testreduce()
        testhelm()
        testfromfile()
        testbadsetup()
        testuserbuffer()
