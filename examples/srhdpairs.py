#!/usr/bin/env python

import optparse
import numpy as np
import matplotlib.pyplot as plt
import h5py
import cowpy
from cowpy import srhdpack


def runstats(fname, opts):
    V = cowpy.fromfile(fname, "prim", members=["vx","vy","vz"], vec3d=True,
                       guard=2, downsample=opts.downsample)
    histpro, histlab = srhdpack.relative_lorentz_pairs(V, 1e6)
    if V.domain.cart_rank == 10:
        import matplotlib.pyplot as plt
        plt.plot(histpro.binloc, histpro.binval, label="propper")
        plt.plot(histlab.binloc, histlab.binval, label="lab frame")
        plt.legend()
        plt.show()

    histpro.name = "gamma-rel-drprop-hist"
    histlab.name = "gamma-rel-drlab-hist"

    if opts.output:
        histpro.dump(opts.output)
        histlab.dump(opts.output)


if __name__ == "__main__":
    parser = optparse.OptionParser()
    parser.add_option("-d", "--downsample", type=int, default=0)
    parser.add_option("-o", "--output", type=str, default="")
    opts, args = parser.parse_args()
    for arg in args:
        runstats(arg, opts)
