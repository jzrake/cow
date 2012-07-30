
import sys
import numpy as np
import cowpy

def test1():
    hist = cowpy.Histogram1d(0, 1)
    for i in range(10):
        hist.add_sample(np.random.rand())
    hist.seal()
    print "got %d counts" % hist.counts

if __name__ == "__main__":
    test1()
