# P4: potential flow initialisation (I 11) on a curvilinear grid: river reach with a bump (discharge inflow,
# stage outflow) Cartesian, curvilinear on the same nodes and turned by 30 deg; the initial velocities
# (after the potential solve) per cell index, the turned ones turned back
import sys, numpy as np
from p3lib import *
import p1lib
nx, ny, nz, dx, h = 60, 12, 8, 0.5, 1.5
L, W = nx*dx, ny*dx
control = "C 11 1\nC 12 21\nC 13 21\nC 14 2\nC 15 21\nC 16 3\n\nB 2 %d %d %d\nB 10 0.0 %g 0.0 %g 0.0 %g\n\nM 10 1\nM 20 2\n" % (nx, ny, nz, L, W, h)
base = os.environ.get("XCTRL","").replace(";","\n") + "A 10 5\nA 520 2\nB 60 1\nF 60 %g\nF 62 %g\nI 11 1\nN 47 0.5\nM 10 1\nP 10 0\nW 10 9.0\nW 22 -9.81\nN 45 1\n" % (h, h)
X0, Y0 = 2.7e5, 7.04e6
th = np.radians(30.0)
I, J = np.meshgrid(np.arange(nx+1), np.arange(ny+1), indexing="ij")
s, n_ = I*dx, J*dx
Xr, Yr = X0 + s*np.cos(th) - n_*np.sin(th), Y0 + s*np.sin(th) + n_*np.cos(th)
Z = 0.5*np.exp(-((s-12.0)**2+(n_-3.0)**2)/6.0)
res = {}
Zf = 0*s
for name, grid, g50 in [("cart", (s, n_, Zf), False), ("curv_flat", (s, n_, Zf), True), ("curv_bump", (s, n_, Z), True), ("rot30_bump", (Xr, Yr, Z), True)]:
    ctrl = base + ("G 50 2\n" if g50 else "")
    d = os.path.join(S3, "runs/t10_"+name)
    rc = p1lib.run(d, control, ctrl, grid=grid, steps=1)
    st = read_state(os.path.join(d, "dump", "state_00000000_r0.bin"))
    res[name] = st
    lap = log_lines(d, "lapiter")
    print("%-10s rc=%d %s max |U| %.4f |V| %.4f |W| %.4f" % (name, rc, lap[-1] if lap else "", np.abs(st["U"]).max(), np.abs(st["V"]).max(), np.abs(st["W"]).max()))
r, c = res["curv_flat"], res["cart"]
print("curv_flat vs cart: U %.1e  V %.1e  W %.1e" % (np.abs(r["U"]-c["U"]).max(), np.abs(r["V"]-c["V"]).max(), np.abs(r["W"]-c["W"]).max()))
r, c = res["rot30_bump"], res["curv_bump"]
U = r["U"]*np.cos(th) + r["V"]*np.sin(th)
V = -r["U"]*np.sin(th) + r["V"]*np.cos(th)
print("rot30_bump vs curv_bump: U %.1e  V %.1e  W %.1e" % (np.abs(U-c["U"]).max(), np.abs(V-c["V"]).max(), np.abs(r["W"]-c["W"]).max()))
