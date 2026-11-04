# T4: non-hydrostatic pressure (A 220 1, 2, 3) and Boussinesq (A 220 4) on curvilinear grids:
#  a) dam break in a closed basin: Cartesian kernels = curvilinear kernels on the same nodes = grid
#     rotated by 90 deg (per cell index)
#     (A 220 4: the inversion is solved with a Krylov method on curvilinear grids, by the line
#     solver on the Cartesian grid: curv=cart to the solver tolerance)
#  b) river reach with an island, discharge inflow: logical nodes vs. rotated by 30 deg at UTM
#     coordinates, with the linear reconstruction A 211 0 (all other terms are rotation invariant)
import sys, numpy as np
from p1lib import *

steps = int(sys.argv[1]) if len(sys.argv) > 1 else 60
cases = [("A220_1", "A 220 1\n"), ("A220_2", "A 220 2\n"), ("A220_3", "A 220 3\nA 224 1.159\n"), ("A220_4", "A 220 4\n")]
if len(sys.argv) > 2:
    cases = [c for c in cases if c[0] in sys.argv[2:]]

# a) dam break
n, dx = 60, 0.5
control_a = "B 1 %g\nB 10 0.0 %g 0.0 %g 0.0 0.5\n\nM 10 1\nM 20 2\n" % (dx, n*dx, n*dx)
ctrl_a = "A 10 2\nA 210 3\nA 211 4\nA 240 1\n%s\nF 60 1.0\nF 72 12.0 18.0 12.0 18.0 1.5\n\nN 41 100.0\nN 47 0.3\n\nM 10 1\n\nW 22 -9.81\n"
I, J = np.meshgrid(np.arange(n+1), np.arange(n+1), indexing="ij")
X0, Y0 = I*dx, J*dx
Z = np.zeros_like(X0)
Xr, Yr = 1000.0-(Y0-15.0), 2000.0+(X0-15.0)

# b) river reach
ni, nj = 100, 20
control_b = "C 11 1\nC 12 3\nC 13 3\nC 14 2\n\nB 1 1.0\nB 10 0.0 %g 0.0 %g 0.0 1.0\n\nM 10 1\nM 20 2\n" % (ni, nj)
ctrl_b = "A 10 2\nA 210 3\nA 211 0\nA 240 1\n%s\nB 60 1\n\nF 50 2\nF 60 2.0\nF 62 2.0\n\nG 50 2\n\nN 41 1000.0\nN 47 0.3\n\nM 10 1\n\nW 10 20.0\nW 22 -9.81\n"
Ib, Jb = np.meshgrid(np.arange(ni+1), np.arange(nj+1), indexing="ij")
Xb, Yb = Ib*1.0, Jb*1.0
Zb = 1.2*np.exp(-((Xb-40.0)**2 + (Yb-14.0)**2)/30.0)
a = np.deg2rad(30.0)
Xb30, Yb30 = 5.0e5 + Xb*np.cos(a) - Yb*np.sin(a), 7.0e6 + Xb*np.sin(a) + Yb*np.cos(a)

def d(a, b, f="eta"):
    return np.abs(a[f]-b[f]).max()

for name, key in cases:
    R = {}
    for v, g, grid, c in [("cart", 0, None, ctrl_a % key),
                          ("id", 2, (X0, Y0, Z), ctrl_a % key + "G 50 2\n"),
                          ("rot90", 2, (Xr, Yr, Z), (ctrl_a % key + "G 50 2\n").replace("F 72 12.0 18.0 12.0 18.0 1.5", "F 72 997.0 1003.0 1997.0 2003.0 1.5"))]:
        dd = os.path.join(S, "runs/t4_%s_%s" % (name, v))
        rc = run(dd, control_a, c, grid=grid, steps=steps)
        R[v] = final_state(dd) if rc == 0 else None
        if rc != 0:
            print(name, v, "FAILED rc", rc, log_lines(dd, "!!!")[:3])
    for v, grid in [("bid", (Xb, Yb, Zb)), ("b30", (Xb30, Yb30, Zb))]:
        dd = os.path.join(S, "runs/t4_%s_%s" % (name, v))
        rc = run(dd, control_b, ctrl_b % key, grid=grid, steps=steps)
        R[v] = final_state(dd) if rc == 0 else None
        if rc != 0:
            print(name, v, "FAILED rc", rc, log_lines(dd, "!!!")[:3])
    if all(R[k] is not None for k in R):
        print("%s: max eta %.4f | curv=cart %.2e | rot90=curv %.2e | river rot30 (A 211 0) %.2e (eta range %.3f)  press: %.2e %.2e %.2e"
              % (name, R["cart"]["eta"].max(), d(R["id"], R["cart"]), d(R["rot90"], R["id"]), d(R["b30"], R["bid"]),
                 np.ptp(R["bid"]["eta"]), d(R["id"], R["cart"], "press"), d(R["rot90"], R["cart"], "press"), d(R["b30"], R["bid"], "press")))
