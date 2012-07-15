
import numpy as np
import cowpy

def test():
    domain = cowpy.UnigridDomain([10,10,10], guard=3)
    dfield = cowpy.UnigridDatafield(domain, ["vx", "vy", "vz"])
    np.random.seed(domain.rank)
    sampx = np.random.rand(10,3)
    dfield.value[:,:,:,0] = 2.1
    dfield.value[:,:,:,1] = 2.2
    dfield.value[:,:,:,2] = 2.3
    x, P = dfield.sample(sampx)
    print x, P
    dfield.value[:,:,:,:] = np.random.rand(*dfield.value.shape)
    x, P = dfield.sample(sampx)
    print x, P
    print domain.global_start
    print domain.coordinate([4, 4, 4])

def testio():
    domain = cowpy.UnigridDomain([12,12,16], guard=3)
    dfield = cowpy.UnigridDatafield(domain, ["vx", "vy", "vz"], name="vel")
    dfield.value[:] = 1.0
    dfield.dump("iotest.h5")
    dfield.read("iotest.h5")
    assert dfield.value.all() == 1.0
    print dfield["vx"].shape

if __name__ == "__main__":
    testio()


