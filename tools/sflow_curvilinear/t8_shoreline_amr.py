# T8: lake at rest with an EMERGED island (peak 1.6 m, water level 1.0 m) with mesh refinement: patch boundaries
# through the island and next to its shoreline, two levels, subcycling, NH, Boussinesq, regridding along the
# shoreline (G 22), two ranks.  The coarse grid does not resolve the shoreline: a coarse cell can be dry while
# some of its children lie below the water.  The velocities must stay at the level of the single-level
# wetting and drying (5e-7 at dx 0.5, 1.3e-4 at dx 0.125: the film A 244 at the shoreline cells)
#
# usage: python3 t8_shoreline_amr.py [steps] [case ...]   (REEF3D_BIN, DIVEMESH_BIN; runs go to ./runs/t8_*)
import sys, os, numpy as np
from p1lib import *

steps = int(sys.argv[1]) if len(sys.argv) > 1 else 400
n, dx = 60, 0.5
L = n*dx
I, J = np.meshgrid(np.arange(n+1), np.arange(n+1), indexing="ij")
X, Y = I*dx, J*dx
Z = 1.6*np.exp(-((X-17.0)**2+(Y-17.5)**2)/12.0)

base = ("A 10 2\nA 210 3\nA 211 4\nA 240 1\nA 220 0\nF 60 1.0\nN 41 100.0\nN 47 0.3\nM 10 %d\nW 22 -9.81\n"
        "G 50 2\nG 3 2\nG 4 8\n")
box = "G 1 1\nG 2 0\nG 10 %s\n"
cases = {
    "single":   ("", 1),                                               # no refinement: the reference
    "cut":      (box % "6.0 17.0 8.0 22.0", 1),                        # patch boundary through the island
    "cut_y":    (box % "6.0 28.0 17.5 28.0", 1),
    "edge":     (box % "6.0 20.0 8.0 22.0", 1),                        # patch boundary 0.6 m off the shoreline
    "all":      (box % "1.0 29.0 1.0 29.0", 1),                        # the whole island refined
    "sub":      (box % "6.0 17.0 8.0 22.0" + "G 7 1\n", 1),
    "lev2":     (box % "6.0 17.0 8.0 22.0" + "G 1 2\n", 1),
    "nh":       (box % "6.0 17.0 8.0 22.0" + "A 220 1\n", 1),
    "bous":     (box % "6.0 17.0 8.0 22.0" + "A 220 4\n", 1),
    "shore":    ("G 1 1\nG 2 5\nG 22 1\n", 1),                          # regridded along the shoreline
    "shore2_sub": ("G 1 2\nG 2 5\nG 22 1\nG 7 1\n", 1),
    "mpi2":     (box % "6.0 17.0 8.0 22.0", 2),
    "shore_mpi2": ("G 1 1\nG 2 5\nG 22 1\n", 2),
}
sel = sys.argv[2:] or list(cases)
print("%-12s %-4s %-10s %-10s %s" % ("case", "rc", "max|U|", "max|V|", "volume change"))
for k in sel:
    extra, np_ = cases[k]
    control = "B 1 %g\nB 10 0.0 %g 0.0 %g 0.0 0.5\n\nM 10 %d\nM 20 2\n" % (dx, L, L, np_)
    d = os.path.join(S, "runs", "t8_" + k)
    rc = run(d, control, (base % np_) + extra, grid=(X, Y, Z), np_=np_, steps=steps)
    um = log_lines(d, r"^umax")[-1].split()[1]
    vm = log_lines(d, r"^vmax")[-1].split()[1]
    log = os.path.join(d, "REEF3D_SFLOW_AMR", "REEF3D_SFLOW_AMR_log.dat")
    dv = open(log).read().splitlines()[-1].split()[-1] if os.path.exists(log) else "-"
    print("%-12s %-4d %-10s %-10s %s" % (k, rc, um, vm, dv))
