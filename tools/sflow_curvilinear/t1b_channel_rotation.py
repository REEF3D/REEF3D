# T1b: river channel with a bed bump (an island near the left bank), discharge inflow and stage
# outflow, on the logical grid itself and rotated by an arbitrary angle (30 deg) and shifted:
# eta and WL per cell index must agree to round-off. A flat-bed Cartesian run checks the
# curvilinear path against the Cartesian kernels.
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
A 211 4
A 220 0
A 240 1

B 60 1

F 50 2
F 60 2.0
F 62 2.0

N 41 1000.0
N 47 0.3

M 10 1

W 10 20.0
W 22 -9.81
"""
steps = int(sys.argv[1]) if len(sys.argv) > 1 else 400

I, J = np.meshgrid(np.arange(ni+1), np.arange(nj+1), indexing="ij")
X0, Y0 = I*dx, J*dx
Zflat = np.zeros_like(X0)
Zbump = 2.3*np.exp(-((X0-40.0)**2 + (Y0-14.0)**2)/30.0)

a = np.deg2rad(30.0)
Xr = 5.0e5 + X0*np.cos(a) - Y0*np.sin(a)
Yr = 7.0e6 + X0*np.sin(a) + Y0*np.cos(a)

res = {}
runs = [("cart_flat", 0, None), ("id_flat", 2, (X0, Y0, Zflat)),
        ("id_bump", 2, (X0, Y0, Zbump)), ("rot_bump", 2, (Xr, Yr, Zbump))]
for name, g, grid in runs:
    c = ctrl + ("G 50 %d\n" % g if g else "")
    d = os.path.join(S, "runs/t1b_"+name)
    rc = run(d, control, c, grid=grid, steps=steps)
    st = final_state(d)
    res[name] = st
    q = log_lines(d, r"^(Inflow|Outflow):")[-2:]
    print("%-10s rc=%d count=%d t=%.4f  eta %.5f..%.5f  %s" % (name, rc, st["count"], st["simtime"], st["eta"].min(), st["eta"].max(), " | ".join(q)))

def cmp(a, b):
    for f in ["eta", "WL"]:
        d = np.abs(res[a][f]-res[b][f]).max()
        print("  %s vs %s: max |d %s| = %.3e" % (a, b, f, d))
    print("  simtime diff %.3e" % (res[a]["simtime"]-res[b]["simtime"]))
cmp("id_flat", "cart_flat")
cmp("rot_bump", "id_bump")
