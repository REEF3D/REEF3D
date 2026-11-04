# T1: dam break in a closed square basin: Cartesian, curvilinear on the same nodes, and
# curvilinear on the grid rotated by 90 deg and shifted; eta per cell index must agree to round-off
import sys, numpy as np
from p1lib import *

n, dx = 120, 0.25
control = """B 1 %g
B 10 0.0 %g 0.0 %g 0.0 0.25

M 10 1
M 20 2
""" % (dx, n*dx, n*dx)

ctrl = """A 10 2
A 210 3
A 211 4
A 220 0
A 240 1

F 60 1.0
F 72 12.0 18.0 12.0 18.0 2.0

N 41 3.0
N 47 0.3

M 10 1

W 22 -9.81
"""
steps = int(sys.argv[1]) if len(sys.argv) > 1 else 80

I, J = np.meshgrid(np.arange(n+1), np.arange(n+1), indexing="ij")
X0, Y0 = I*dx, J*dx
Z = np.zeros_like(X0)

# exact rotation by 90 deg about (15,15), shifted: the dam box maps to a box
Xr = 1000.0 - (Y0-15.0)
Yr = 2000.0 + (X0-15.0)

res = {}
for name, g, grid in [("cart", 0, None), ("curv_id", 2, (X0, Y0, Z)), ("curv_rot", 2, (Xr, Yr, Z))]:
    c = ctrl + ("G 50 %d\n" % g if g else "")
    if name == "curv_rot":
        c = c.replace("F 72 12.0 18.0 12.0 18.0 2.0", "F 72 997.0 1003.0 1997.0 2003.0 2.0")
    rc = run(os.path.join(S, "runs/t1_"+name), control, c, grid=grid, steps=steps)
    st = final_state(os.path.join(S, "runs/t1_"+name))
    res[name] = st
    print("%-9s rc=%d count=%d t=%.6f  max eta %.6f" % (name, rc, st["count"], st["simtime"], st["eta"].max()))

for name in ["curv_id", "curv_rot"]:
    for f in ["eta", "WL"]:
        d = np.abs(res[name][f]-res["cart"][f]).max()
        print("%-9s vs cart: max |d %s| = %.3e  (max %s %.3f)" % (name, f, d, f, np.abs(res["cart"][f]).max()))
    print("%-9s simtime diff %.3e" % (name, res[name]["simtime"]-res["cart"]["simtime"]))
