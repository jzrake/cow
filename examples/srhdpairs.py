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
    histpro, histlab = srhdpack.relative_lorentz_pairs(V, opts.pairs)
    histpro.name = "gamma-rel-drprop-hist"
    histlab.name = "gamma-rel-drlab-hist"

    if opts.output:
        histpro.dump(opts.output)
        histlab.dump(opts.output)

    if opts.plot and V.domain.cart_rank == 0:
        import matplotlib.pyplot as plt
        plt.plot(histpro.binloc, histpro.binval, label=histpro.name)
        plt.plot(histlab.binloc, histlab.binval, label=histlab.name)
        plt.legend()
        plt.show()


if __name__ == "__main__":
    parser = optparse.OptionParser()
    parser.add_option("-p", "--pairs", type=float, default=1e6)
    parser.add_option("-d", "--downsample", type=int, default=0)
    parser.add_option("-o", "--output", type=str, default="")
    parser.add_option("--plot", action="store_true", default=False)
    opts, args = parser.parse_args()
    for arg in args:
        runstats(arg, opts)
