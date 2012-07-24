#!/usr/bin/env python

import optparse
import numpy as np
import cowpy
from cowpy import relative_lorentz_pairs as pairs


def cversion(fname, opts):
    V = cowpy.fromfile(fname, "prim", members=["vx","vy","vz"], vec3d=True,
                       guard=2, downsample=opts.downsample)
    histpro, histlab = pairs(V, opts.pairs, bins=opts.bins)
    histpro.name = "gamma-rel-drprop-hist"
    histlab.name = "gamma-rel-drlab-hist"

    if opts.output:
        histpro.dump(opts.output)
        histlab.dump(opts.output)

    if opts.plot and V.domain.cart_rank == 0:
        import matplotlib.pyplot as plt
        plt.plot(histpro.binloc, histpro.binval, label=histpro.name)
        plt.plot(histlab.binloc, histlab.binval, label=histlab.name)
        plt.legend(loc='best')
        plt.ylim(1.0, 2.8)
        plt.show()


def gammarel(fname, opts):
    V = cowpy.fromfile(fname, "prim", members=["vx","vy","vz"], vec3d=True,
                       guard=2, downsample=opts.downsample)
    U = cowpy.DataField(V.domain, members=["u0", "u1", "u2", "u3"])
    np.random.seed(V.domain.cart_rank)
    nsamp = int(opts.pairs * 2)
    npair = int(opts.pairs)
    sampx = np.random.rand(nsamp, 3)

    vx, vy, vz = V.value[...,0], V.value[...,1], V.value[...,2]
    U.value[...,0] = 1.0 / (1.0 - (vx**2 + vy**2 + vz**2))**0.5
    U.value[...,1] = vx * U.value[...,0]
    U.value[...,2] = vy * U.value[...,0]
    U.value[...,3] = vz * U.value[...,0]
    U.sync_guard()

    x, P = U.sample_global(sampx)
    histP = cowpy.Histogram1d(0.0, 1.5, bins=opts.bins, binmode="average")
    for n in range(npair):
        i0 = np.random.randint(0, len(P))
        i1 = np.random.randint(0, len(P))
        dx = np.dot(x[i1] - x[i0], x[i1] - x[i0])**0.5
        dg = P[i0][0]*P[i1][0] - np.dot(P[i0][1:], P[i1][1:])
        histP.add_sample(dx, weight=dg)
    histP.seal()
    histP.name = "gamma-rel-drlab-hist"

    if opts.output:
        histP.dump(opts.output)

    if opts.plot and V.domain.cart_rank == 0:
        import matplotlib.pyplot as plt
        plt.plot(histP.binloc, histP.binval, label="P")
        plt.legend(loc='best')
        plt.ylim(1.0, 2.8)
        plt.show()


def vrel(fname, opts):
    V = cowpy.fromfile(fname, "prim", members=["vx","vy","vz"], vec3d=True,
                       guard=2, downsample=opts.downsample)
    np.random.seed(V.domain.cart_rank)
    nsamp = int(opts.pairs * 2)
    npair = int(opts.pairs)
    sampx = np.random.rand(nsamp, 3)
    x, P = V.sample_global(sampx)
    histP = cowpy.Histogram1d(0.0, 1.5, bins=opts.bins, binmode="average")

    for n in range(npair):
        i0 = np.random.randint(0, len(P))
        i1 = np.random.randint(0, len(P))
        dx = np.dot(x[i1] - x[i0], x[i1] - x[i0])**0.5
        dP = np.dot(P[i1] - P[i0], P[i1] - P[i0])**0.5
        histP.add_sample(dx, weight=dP)
    histP.seal()

    if V.domain.cart_rank == 0:
        import matplotlib.pyplot as plt
        plt.plot(histP.binloc, histP.binval, label="P")
        plt.legend(loc='best')
        plt.show()


if __name__ == "__main__":
    parser = optparse.OptionParser()
    parser.add_option("-p", "--pairs", type=float, default=1e6)
    parser.add_option("-d", "--downsample", type=int, default=0)
    parser.add_option("-o", "--output", type=str, default="")
    parser.add_option("--plot", action="store_true", default=False)
    parser.add_option("--bins", type=int, default=72)
    opts, args = parser.parse_args()
    for arg in args:
        cversion(arg, opts)
        #vrel(arg, opts)
        #gammarel(arg, opts)
