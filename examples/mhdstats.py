#!/usr/bin/env python

import optparse
import numpy as np
import matplotlib.pyplot as plt
import h5py
import cowpy

def maghist(field):
    d = field.domain
    min, max, sum = field.reduce_magnitude()
    hist = cowpy.Histogram1d(min, max, bins=100, binmode="counts", domain=d)
    N = len(field.members)
    if N == 3:
        mag = np.sqrt(field.interior[...,0]**2 +
                      field.interior[...,1]**2 + field.interior[...,2]**2)
    elif N == 1:
        mag = field.interior[...,0]
    else:
        raise ValueError()
    for a in mag.flatten():
        hist.add_sample(a)
    hist.seal()
    hist.name = field.name + "-hist"
    return hist

def runstats(fname, opts):
    V = cowpy.fromfile(fname, "prim", members=["vx","vy","vz"], vec3d=True,
                       guard=2, downsample=opts.downsample)
    B = cowpy.fromfile(fname, "prim", members=["Bx","By","Bz"], vec3d=True,
                       guard=2, downsample=opts.downsample)
    V.name = "V"
    B.name = "B"
    f = { }
    f["divV"] = V.divergence()
    f["divB"] = B.divergence()
    f["curlV"] = V.curl()
    f["curlB"] = B.curl()
    f["vdotB"] = cowpy.dot_product(V, B)
    f["vcrossB"] = cowpy.cross_product(V, B)
    f["curlBdotvcrossB"] = cowpy.cross_product(f["curlB"], f["vcrossB"])
    f["curlBdotB"] = cowpy.dot_product(f["curlB"], B)
    f["vcrossBcrossB"] = cowpy.cross_product(f["vcrossB"], B)
    f["divvcrossBcrossB"] = f["vcrossBcrossB"].divergence()

    for thing in f:
        hist = maghist(f[thing])
        hist.name = thing + "-hist"
        hist.dump(opts.output or fname)

if __name__ == "__main__":
    parser = optparse.OptionParser()
    parser.add_option("-d", "--downsample", type=int, default=0)
    parser.add_option("-o", "--output", type=str, default="")
    opts, args = parser.parse_args()
    for arg in args:
        runstats(arg, opts)
