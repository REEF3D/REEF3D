# T9: horizontal diffusion on the curvilinear grid (A 212 1, sflow_ediff::curv_stress) with a constant viscosity
# (W 2) and the algebraic eddy viscosity (A 260 3 mixing length, A 260 4 parabolic): the river channel of T1b with an
# emerged island near the left bank (a wake with shear layers), on the logical grid itself and rotated by 30 deg and
# shifted (A 211 0: the WENO weights act on the Cartesian components, see t1b), and on a skewed grid.  The stress is a tensor of the Cartesian velocity gradients: eta, WL per cell index of
# the identity and the rotated grid agree to round-off; the diffusion changes the wake (against no diffusion).
# usage: python3 t9_diffusion_rotation.py [steps]   (REEF3D_BIN, DIVEMESH_BIN; runs go to ./runs/t9_*)
import sys, numpy as np
from p1lib import *

ni, nj, dx = 100, 20, 1.0
control = """C 11 1
C 12 3
C 13 3
C 14 2

B 1 %g
B 10 0.0 %g 0.0 %g 0.0 1.0

M 10 1
M 20 2
""" % (dx, ni*dx, nj*dx)

ctrl = """A 10 2
A 210 3
A 211 0
A 220 0
A 240 1
B 50 0.01
B 60 1
F 50 2
F 60 2.0
F 62 2.0
G 50 2
N 41 1000.0
N 47 0.3
M 10 1
W 10 20.0
W 22 -9.81
"""
steps = int(sys.argv[1]) if len(sys.argv) > 1 else 600

I, J = np.meshgrid(np.arange(ni+1), np.arange(nj+1), indexing="ij")
X0, Y0 = I*dx, J*dx
Z = 2.3*np.exp(-((X0-40.0)**2 + (Y0-14.0)**2)/30.0)
a = np.deg2rad(30.0)
Xr = 5.0e5 + X0*np.cos(a) - Y0*np.sin(a)
Yr = 7.0e6 + X0*np.sin(a) + Y0*np.cos(a)
Xs = X0 + 0.25*Y0 + 0.6*np.sin(2*np.pi*X0/50.0)*np.sin(np.pi*Y0/20.0)   # skewed, wavy
Ys = Y0 + 0.4*np.sin(2*np.pi*X0/25.0)*np.sin(np.pi*Y0/20.0)

models = {"none": "", "visc": "A 212 1\nW 2 0.05\n", "mixlen": "A 260 3\n", "parab": "A 260 4\nW 2 0.01\n"}
res = {}
for m, extra in models.items():
    grids = [("id", (X0, Y0, Z)), ("rot", (Xr, Yr, Z))] + ([("skew", (Xs, Ys, Z))] if m == "mixlen" else [])
    for g, grid in grids:
        d = os.path.join(S, "runs/t9_%s_%s" % (m, g))
        rc = run(d, control, ctrl + extra, grid=grid, steps=steps)
        st = final_state(d)
        res[(m, g)] = st
        print("%-7s %-5s rc=%d t=%.3f  eta %.5f .. %.5f" % (m, g, rc, st["simtime"], st["eta"].min(), st["eta"].max()))
for m in models:
    a_, b_ = res[(m, "id")], res[(m, "rot")]
    print("%-7s rot vs id: eta %.2e  WL %.2e  (simtime %.1e)" % (m, np.abs(a_["eta"]-b_["eta"]).max(),
          np.abs(a_["WL"]-b_["WL"]).max(), a_["simtime"]-b_["simtime"]))
for m in models:
    if m != "none":
        print("%-7s vs none (id): eta %.2e" % (m, np.abs(res[(m, "id")]["eta"]-res[("none", "id")]["eta"]).max()))
