# Weber, Schumate & Mawer (2001) 90 deg open-channel junction (Shumate 1998, IIHR flume) as a three-block grid:
#   block 0  the main channel upstream of the junction (native inflow Qm on its -i face, B 60 1, W 10)
#   block 1  the junction square and the main channel downstream of it (water level H0 on its +i face, G 54 1 2)
#   block 2  the branch channel, i towards the junction (inflow Qb on its -i face, G 54 2 1)
# The flume: main channel 21.946 m (72 ft), branch 3.658 m (12 ft) joining 5.486 m (18 ft) downstream of the main
# inlet, both 0.914 m (W) wide, horizontal bed, Qt = 0.170 m3/s, tailwater 0.296 m.  Layout: the main channel along
# +x over y 0 .. W, the branch from +y (y W .. W + Lb, x Lu .. Lu + W).  The coordinates of Huang et al. (2002):
# x* = (Lu - x)/W (positive upstream, 0 at the upstream corner of the branch), y* = (W - y)/W (0 at the bank of
# the branch), so the branch lies at -1 < x* < 0, y* < 0.
# Writes blocks.json (LAGOON's block topology), blocks.dat, control.txt and ctrl.txt into the case folder.
# usage: python3 weber_case.py casedir [--qstar 0.75] [--n 16] [--solver sflow|nhflow] [--np 2] [--tend 200]
#        [--out 20] [--nz 10] [--ks 0.001] [--xctrl "A 560 1;..."]  (n: cells across W)
import argparse, json, os, numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("case")
ap.add_argument("--qstar", type=float, default=0.75, help="Qm / Qt (RUN2 0.25, RUN5 0.75)")
ap.add_argument("--qt", type=float, default=0.170)
ap.add_argument("--h0", type=float, default=0.296, help="water level at the outlet (the tailgate)")
ap.add_argument("--n", type=int, default=16, help="cells across a channel")
ap.add_argument("--solver", default="sflow")
ap.add_argument("--np", type=int, default=2)
ap.add_argument("--tend", type=float, default=200.0)
ap.add_argument("--out", type=float, default=20.0)
ap.add_argument("--nz", type=int, default=10)
ap.add_argument("--ks", type=float, default=0.001, help="bed roughness ks (B 50)")
ap.add_argument("--xctrl", default="")
ap.add_argument("--layout", default="a", help="a: native inflow on block 0, b: native outflow on block 0")
a = ap.parse_args()
os.makedirs(a.case, exist_ok=True)

FT = 0.3048
W, Lm, Lb, Lu = 3*FT, 72*FT, 12*FT, 18*FT
Ld = Lm - Lu - W
n = a.n
dx = W/n
nu, nd, nb = round(Lu/dx), round(Ld/dx), round(Lb/dx)
Qm, Qb = a.qstar*a.qt, (1.0 - a.qstar)*a.qt

# the nodes of the blocks (i, j), straight Cartesian lines
def box(x0, y0, ex, ey, ni, nj):
    I, J = np.meshgrid(np.arange(ni+1), np.arange(nj+1), indexing="ij")
    X = x0 + ex[0]*I*dx + ey[0]*J*dx
    Y = y0 + ex[1]*I*dx + ey[1]*J*dx
    return X, Y, np.zeros_like(X)
V = {"A": (0.0, 0.0), "B": (0.0, W), "C": (Lu, 0.0), "D": (Lu, W), "E": (Lu+W, W), "F": (Lm, 0.0), "G": (Lm, W),
     "H": (Lu, W+Lb), "I": (Lu+W, W+Lb), "K": (Lu+W, 0.0)}
def edge(s, e, c):
    return {"start": s, "end": e, "cells": c, "curve": "straight"}
branch = {"name": "branch", "faces": {"i0": ["branch_inlet"], "i1": ["mouth"], "j0": ["branch_up"], "j1": ["branch_down"]},
          "bc": {"i0": {"type": "inflow", "value": round(Qb, 6)}}, "smoothing": 0}
edges = {"inlet": edge("A", "B", n), "up_left": edge("B", "D", nu), "mouth": edge("D", "E", n),
         "down_left": edge("E", "G", nd), "outlet": edge("F", "G", n),
         "branch_inlet": edge("H", "I", n), "branch_up": edge("H", "D", nb), "branch_down": edge("I", "E", nb)}
if a.layout == "a":
    # block 0 the main channel upstream of the junction (native inflow), block 1 the junction and the main channel
    # downstream (water level on its +i face)
    edges.update({"up_right": edge("A", "C", nu), "cut": edge("C", "D", n), "junction_right": edge("C", "F", n+nd)})
    blocks = [
        {"name": "main upstream", "faces": {"i0": ["inlet"], "i1": ["cut"], "j0": ["up_right"], "j1": ["up_left"]},
         "bc": {"i0": {"type": "inflow", "value": round(Qm, 6)}}, "smoothing": 0},
        {"name": "junction and main downstream", "faces": {"i0": ["cut"], "i1": ["outlet"], "j0": ["junction_right"],
                                                           "j1": ["mouth", "down_left"]},
         "bc": {"i1": {"type": "outflow", "value": a.h0}}, "smoothing": 0},
        branch]
    grids = [box(0.0, 0.0, (1, 0), (0, 1), nu, n),           # i +x, j +y
             box(Lu, 0.0, (1, 0), (0, 1), n+nd, n),
             box(Lu, W+Lb, (0, -1), (1, 0), nb, n)]          # i -y, j +x
else:
    # block 0 the main channel downstream of the junction (native outflow), block 1 the main channel upstream and the
    # junction (inflow on its -i face)
    edges.update({"up_right": edge("A", "K", nu+n), "cut": edge("K", "E", n), "down_right": edge("K", "F", nd)})
    blocks = [
        {"name": "main downstream", "faces": {"i0": ["cut"], "i1": ["outlet"], "j0": ["down_right"], "j1": ["down_left"]},
         "bc": {"i1": {"type": "outflow", "value": a.h0}}, "smoothing": 0},
        {"name": "main upstream and junction", "faces": {"i0": ["inlet"], "i1": ["cut"], "j0": ["up_right"],
                                                         "j1": ["up_left", "mouth"]},
         "bc": {"i0": {"type": "inflow", "value": round(Qm, 6)}}, "smoothing": 0},
        branch]
    grids = [box(Lu+W, 0.0, (1, 0), (0, 1), nd, n),
             box(0.0, 0.0, (1, 0), (0, 1), nu+n, n),
             box(Lu, W+Lb, (0, -1), (1, 0), nb, n)]
used = {e[k] for e in edges.values() for k in ("start", "end")}
topo = {"version": 1, "vertices": {k: [round(x, 6), round(y, 6)] for k, (x, y) in V.items() if k in used}, "edges": edges,
        "blocks": blocks, "bed": "flat", "z": 0.0,
        "note": "Weber et al. (2001) 90 deg junction, q* = %g (Qm %.4f, Qb %.4f m3/s), W/dx = %d, layout %s" % (a.qstar, Qm, Qb, n, a.layout)}
json.dump(topo, open(os.path.join(a.case, "blocks.json"), "w"), indent=1)

with open(os.path.join(a.case, "blocks.dat"), "w") as f:
    f.write("%d\n" % len(grids))
    for X, Y, Z in grids:
        f.write("%d %d\n" % (X.shape[0]-1, X.shape[1]-1))
        for j in range(X.shape[1]):
            for i in range(X.shape[0]):
                f.write("%.12g %.12g %.6g\n" % (X[i, j], Y[i, j], Z[i, j]))

H = a.nz*dx if a.solver == "nhflow" else 0.6
c11, c14 = (1, 21) if a.layout == "a" else (21, 2)
open(os.path.join(a.case, "control.txt"), "w").write(
    "// Weber et al. (2001) 90 deg junction: the three blocks of blocks.json\nC 11 %d\nC 12 21\nC 13 21\nC 14 %d\nC 15 21\n"
    "C 16 3\n\nB 1 %g\nB 10 %g %g 0.0 %g 0.0 %g\n\nM 10 %d\n\nR 20 1\n" % (c11, c14, dx, 0.0, Lm, W+Lb, H, a.np))
if a.solver == "sflow":
    # A 218 1: the bed friction of SFLOW (Manning, n = ks^(1/6)/20); B 50 the ks of the same Chezy at h0 as the log law of
    # NHFLOW with the Nikuradse roughness --ks (C = 18 log10(12 h0/ks))
    flow = "A 10 2\nA 210 3\nA 211 4\nA 218 1\nA 220 0\nA 240 1\nP 181 -1\nP 182 %g\n" % a.out
else:
    # A 519 1: the bed friction of the bottom layer (without it the NHFLOW bed is frictionless, B 50 alone acts only
    # on the turbulence wall functions); B 61 2: the log-law profile at the discharge inflows (native and G 54)
    # A 567 3: the free-surface value of epsilon/omega (Celik & Rodi, y' = 0.07 h): without it k-epsilon and k-omega keep
    # their eddy viscosity up to the surface (zero gradient) and mix the shear layer of the separation over the depth
    flow = "A 10 5\nA 519 1\nA 520 0\nA 567 3\nB 61 2\nP 10 1\nP 30 %g\n" % a.out
import math
ks_b50 = a.ks if a.solver == "nhflow" else (20.0*a.h0**(1.0/6.0)/(18.0*math.log10(12.0*a.h0/a.ks)))**6
g54 = ("G 54 1 2 %g\nG 54 2 1 %.6f\n" % (a.h0, Qb) if a.layout == "a" else "G 54 1 1 %.6f\nG 54 2 1 %.6f\n" % (Qm, Qb))
ctrl = ("// Weber et al. (2001) 90 deg junction, q* = %g, layout %s\n" % (a.qstar, a.layout) + flow +
        "B 50 %g\nB 60 1\nW 10 %.6f\nF 50 2\nF 60 %g\nF 62 %g\nG 50 1\n%sN 41 %g\nN 47 0.3\n"
        "M 10 %d\nW 22 -9.81\n" % (ks_b50, Qm, a.h0, a.h0, g54, a.tend, a.np) + a.xctrl.replace(";", "\n") + "\n")
open(os.path.join(a.case, "ctrl.txt"), "w").write(ctrl)
print("W %.4f dx %.5f cells: upstream %d x %d, junction + downstream %d x %d, branch %d x %d (%d); Qm %.4f Qb %.4f"
      % (W, dx, nu, n, n+nd, n, nb, n, (nu+n+nd+nb)*n, Qm, Qb))
