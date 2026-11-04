# The Weber junction of weber_case.py on ONE Cartesian grid with native SFLOW, the reference for the blocks and for
# the turbulence models: the box of the flume (x 0 .. 72 ft, y 0 .. W + 12 ft) with DIVEMesh solids (S 10) beside the
# branch, the main inflow native (C 11 1, B 60 1, W 10 Qm), the branch inflow a patch on the +y face (B 440 face 2,
# B 411 -Qb: the discharge of a patch counts along +y), the water level native at +x (C 14 2, F 62).  The same grid lines as the blocks (dx = W/n).
# usage: python3 weber_cart.py casedir [--qstar 0.75] [--n 16] [--np 2] [--tend 600] [--turb 0] [--xctrl "W 2 1e-3"]
#        (turb: A 260; 0 none, 1 depth-averaged k-eps, 2 k-omega, 3 Prandtl mixing length, 4 parabolic)
import argparse, json, os

ap = argparse.ArgumentParser()
ap.add_argument("case")
ap.add_argument("--qstar", type=float, default=0.75)
ap.add_argument("--qt", type=float, default=0.170)
ap.add_argument("--h0", type=float, default=0.296)
ap.add_argument("--n", type=int, default=16)
ap.add_argument("--np", type=int, default=2)
ap.add_argument("--tend", type=float, default=600.0)
ap.add_argument("--out", type=float, default=100.0)
ap.add_argument("--ks", type=float, default=0.001)
ap.add_argument("--turb", type=int, default=0)
ap.add_argument("--xctrl", default="")
a = ap.parse_args()
os.makedirs(a.case, exist_ok=True)
FT = 0.3048
W, Lm, Lb, Lu = 3*FT, 72*FT, 12*FT, 18*FT
dx = W/a.n
Qm, Qb = a.qstar*a.qt, (1.0 - a.qstar)*a.qt
e = 0.25*dx
open(os.path.join(a.case, "control.txt"), "w").write(
    "// Weber et al. (2001) 90 deg junction, one Cartesian grid (weber_cart.py)\nC 11 1\nC 12 21\nC 13 21\nC 14 2\n"
    "C 15 21\nC 16 3\n\nB 1 %g\nB 10 0.0 %g 0.0 %g 0.0 0.6\n\nS 10 %g %g %g %g -1.0 2.0\nS 10 %g %g %g %g -1.0 2.0\n\n"
    "M 10 %d\n" % (dx, Lm, W+Lb, -1.0, Lu, W, W+Lb+1.0, Lu+W, Lm+1.0, W, W+Lb+1.0, a.np))
turb = "A 260 %d\n" % a.turb if a.turb else ""
# A 218 1: the bed friction (Manning, n = ks^(1/6)/20); B 50 the ks of the same Chezy at h0 as the log law with the
# Nikuradse roughness --ks (C = 18 log10(12 h0/ks)), as weber_case.py
import math
ks_b50 = (20.0*a.h0**(1.0/6.0)/(18.0*math.log10(12.0*a.h0/a.ks)))**6
ctrl = ("// Weber et al. (2001) 90 deg junction, q* = %g, one Cartesian grid\nA 10 2\nA 210 3\nA 211 4\nA 218 1\nA 220 0\nA 240 1\n"
        "%sB 50 %g\nB 60 1\nW 10 %.6f\nB 440 1 2 %g %g %g %g\nB 411 1 %.6f\nF 50 2\nF 60 %g\nF 62 %g\nN 41 %g\nN 47 0.3\n"
        "M 10 %d\nP 181 -1\nP 182 %g\nW 22 -9.81\n" % (a.qstar, turb, ks_b50, Qm, Lu-e, Lu+W+e, W+Lb-e, W+Lb+e, -Qb, a.h0,
                                                     a.h0, a.tend, a.np, a.out) + a.xctrl.replace(";", "\n") + "\n")
open(os.path.join(a.case, "ctrl.txt"), "w").write(ctrl)
json.dump({"weber_cart": 1, "qstar": a.qstar, "n": a.n, "dx": dx,
           "note": "Weber et al. (2001) 90 deg junction, q* = %g, one Cartesian grid" % a.qstar},
          open(os.path.join(a.case, "case.json"), "w"), indent=1)
