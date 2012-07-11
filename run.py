

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


if __name__ == "__main__":
    test()


