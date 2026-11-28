#!/usr/bin/env python3
# staged_release.py: the release times of a staged (retrogressive) landslide for SFLOW F 74 (sflow-release.dat),
# from a case made by slide_terrain.py (reference.npz: the release mask and the pre / post terrain).
#
#   python3 staged_release.py <case dir> --blocks 4 --dt 60 [--dir 90] [--t0 0] [--equal volume|length] [--png]
#
# - the release area is cut into --blocks bands across the direction --dir (degrees from east, counter-clockwise;
#   90 = the retrogression runs to the north: the first block is the southern toe, the last the backscarp)
# - --equal volume (default): the bands hold equal volumes (pre - post); --equal length: equal widths along --dir
# - block k (0 = toe) is released at t0 + k*dt; block 0 at t0 = 0 is in the flow from the start
# - sflow-release.dat has the grid of sflow-level.dat (nx ny x0 y0 dx dy, then nx*ny values, i fastest); outside the
#   release: -9999 (released at the start; there is no material of F 73 there anyway)
# - --png: a map of the blocks (staged_release.png)

import argparse
import json
import os

import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case")
    ap.add_argument("--blocks", type=int, default=4)
    ap.add_argument("--dt", type=float, required=True)
    ap.add_argument("--t0", type=float, default=0.0)
    ap.add_argument("--dir", type=float, default=90.0)
    ap.add_argument("--equal", choices=("volume", "length"), default="volume")
    ap.add_argument("--png", action="store_true")
    a = ap.parse_args()

    r = np.load(os.path.join(a.case, "reference.npz"))
    x, y, rel = r["x"], r["y"], r["release"]
    thick = np.where(rel, r["pre"] - r["post"], 0.0)
    dx = float(r["dx"])
    X, Y = np.meshgrid(x, y)
    th = np.radians(a.dir)
    s = X*np.cos(th) + Y*np.sin(th)                 # the coordinate along the retrogression

    sr, vr = s[rel], thick[rel]
    order = np.argsort(sr)
    if a.equal == "volume":
        cv = np.cumsum(vr[order])/vr.sum()
        edges = [sr[order][np.searchsorted(cv, k/a.blocks)] for k in range(1, a.blocks)]
    else:
        edges = list(np.linspace(sr.min(), sr.max(), a.blocks + 1)[1:-1])
    block = np.searchsorted(np.array(edges), s, side="right")

    t = np.where(rel, a.t0 + block*a.dt, -9999.0)
    ny, nx = t.shape
    with open(os.path.join(a.case, "sflow-release.dat"), "w") as f:
        f.write("%d %d %.3f %.3f %.3f %.3f\n" % (nx, ny, x[0], y[0], dx, dx))
        for j in range(ny):
            f.write(" ".join("%.1f" % v for v in t[j]) + "\n")

    info = []
    for k in range(a.blocks):
        m = rel & (block == k)
        info.append(dict(block=k, time=a.t0 + k*a.dt, cells=int(m.sum()), volume_m3=float(thick[m].sum()*dx*dx)))
        print("block %d: t = %6.1f s, %5d cells, %.3f Mm3" % (k, a.t0 + k*a.dt, m.sum(), thick[m].sum()*dx*dx/1e6))
    with open(os.path.join(a.case, "staged_release.json"), "w") as f:
        json.dump(dict(blocks=a.blocks, dt=a.dt, t0=a.t0, dir=a.dir, equal=a.equal, info=info), f, indent=1)

    if a.png:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6, 6))
        jj, ii = np.where(rel)
        sl = (slice(max(jj.min() - 20, 0), jj.max() + 20), slice(max(ii.min() - 20, 0), ii.max() + 20))
        ext = [x[sl[1]][0] - dx/2, x[sl[1]][-1] + dx/2, y[sl[0]][0] - dx/2, y[sl[0]][-1] + dx/2]
        ax.imshow(r["pre"][sl], origin="lower", extent=ext, cmap="gray", alpha=0.6)
        bm = np.ma.masked_where(~rel[sl], block[sl])
        im = ax.imshow(bm, origin="lower", extent=ext, cmap="viridis", vmin=0, vmax=a.blocks - 1)
        cb = fig.colorbar(im, ax=ax, ticks=range(a.blocks), shrink=0.8)
        cb.set_ticklabels(["%d s" % (a.t0 + k*a.dt) for k in range(a.blocks)])
        ax.set_title("staged release: %d blocks, dt %.0f s" % (a.blocks, a.dt))
        ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]")
        fig.savefig(os.path.join(a.case, "staged_release.png"), dpi=120, bbox_inches="tight")


if __name__ == "__main__":
    main()
