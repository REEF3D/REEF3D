#!/usr/bin/env python3
# Architect: Hans Bihs
"""
Stretched rectilinear grid for DIVEMesh (B 101 10 / B 102 10): x-spacing.dat and y-spacing.dat, the node coordinates
of each axis divided by the domain length (DIVEMesh: XN = xmin + value * xmax, so with xmin = 0).

Fine cells in bands, a geometric growth away from them: the spacing at x is

    d(x) = min(dmax, min over the bands (a, b, d_b) of  d_b + (r - 1) dist(x, [a, b]))

so that neighbouring cells differ by about the factor r; the nodes are marched with the spacing at the middle of each
cell and the last cells are stretched to end exactly at the domain length. Useful for the semi-implicit SFLOW
(A 226 1), which has no mesh refinement (G 1): fine cells over a coast, coarse ones offshore.

  python3 stretched_grid.py --lx 460000 --ly 316000 --xband 180000 300000 250 --xband 300000 460000 500 \\
                            --yband 205000 292000 250 --yband 292000 316000 500 [--ratio 1.08] [--dmax 2000] [-o .]

DIVEMesh control.txt:  B 10 0.0 <lx> 0.0 <ly> 0.0 <dz>   B 101 10   B 102 10   (B 1 sets the z cells only)
Bathymetry (geo.dat): one point per cell centre is the safe choice (the centres are printed with --centres).
"""
import argparse
import os

import numpy as np


def axis(L, bands, ratio=1.08, dmax=2000.0):
    """node coordinates 0..L of one axis"""
    def d(x):
        v = dmax
        for a, b, db in bands:
            v = min(v, db + (ratio - 1.0) * max(a - x, x - b, 0.0))
        return v
    x, nodes = 0.0, [0.0]
    while x < L:
        h = d(x)
        x += d(x + 0.5 * h)
        nodes.append(x)
    nodes = np.array(nodes)
    # the overshoot spread over the last cells in proportion to their size
    k = max(2, len(nodes) // 20)
    seg = nodes[-k - 1:] - nodes[-k - 1]
    nodes[-k - 1:] = nodes[-k - 1] + seg * (L - nodes[-k - 1]) / seg[-1]
    return nodes


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--lx", type=float, required=True)
    ap.add_argument("--ly", type=float, default=0.0, help="0: no y-spacing.dat (2D in x only)")
    ap.add_argument("--xband", type=float, nargs=3, action="append", default=[], metavar=("A", "B", "D"))
    ap.add_argument("--yband", type=float, nargs=3, action="append", default=[], metavar=("A", "B", "D"))
    ap.add_argument("--ratio", type=float, default=1.08)
    ap.add_argument("--dmax", type=float, default=2000.0)
    ap.add_argument("--centres", action="store_true", help="also write x-centres.dat, y-centres.dat [m]")
    ap.add_argument("-o", default=".")
    a = ap.parse_args()
    for name, L, bands in (("x", a.lx, a.xband), ("y", a.ly, a.yband)):
        if not L > 0.0:
            continue
        n = axis(L, [tuple(b) for b in bands], a.ratio, a.dmax)
        dx = np.diff(n)
        r = np.max(np.maximum(dx[1:] / dx[:-1], dx[:-1] / dx[1:])) if len(dx) > 1 else 1.0
        np.savetxt(os.path.join(a.o, "%s-spacing.dat" % name), n / L, fmt="%.9f")
        if a.centres:
            np.savetxt(os.path.join(a.o, "%s-centres.dat" % name), 0.5 * (n[1:] + n[:-1]), fmt="%.3f")
        print("%s: %d cells, spacing %.1f .. %.1f m, largest ratio of neighbours %.3f" % (name, len(dx), dx.min(), dx.max(), r))


if __name__ == "__main__":
    main()
