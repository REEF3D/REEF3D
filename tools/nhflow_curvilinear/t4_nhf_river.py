# P3 T4: river reach (B 60 1 discharge inflow, F 62 stage outflow, bed friction): Cartesian, curvilinear
# on the same nodes, and turned by 30 deg at UTM coordinates (per cell index, velocities turned back);
# discharge in and out from the log. VAR: "" | ediff | idiff (horizontal diffusion A 512 1 / 2)
import sys, numpy as np, re
from p3lib import *
import p1lib
steps = int(sys.argv[1]) if len(sys.argv) > 1 else 300
nx, ny, nz, dx, h = 60, 12, 8, 0.5, 1.5
L, W = nx*dx, ny*dx
control = "C 11 1\nC 12 21\nC 13 21\nC 14 2\nC 15 21\nC 16 3\n\nB 2 %d %d %d\nB 10 0.0 %g 0.0 %g 0.0 %g\n\nM 10 1\nM 20 2\n" % (nx, ny, nz, L, W, h)
A520 = int(os.environ.get("A520","2"))
var = os.environ.get("VAR","")
extra = {"": "", "ediff": "A 512 1\nW 2 0.01\n", "idiff": "A 512 2\nW 2 0.01\n", "ke": "A 560 1\n", "kw": "A 560 2\n", "les": "A 560 31\n"}[var]
base = os.environ.get("XCTRL","").replace(";","\n") + "A 10 5\nA 519 1\nA 520 %d\nB 50 0.01\nB 60 1\nF 60 %g\nF 62 %g\nN 47 0.5\nN 48 0\nN 49 0.02\nM 10 1\nP 10 0\nP 12 1\nW 10 9.0\nW 22 -9.81\n%s" % (A520, h, h, extra)
X0, Y0 = 2.7e5, 7.04e6
th = np.radians(float(os.environ.get("ROT","30.0")))
I, J = np.meshgrid(np.arange(nx+1), np.arange(ny+1), indexing="ij")
s, n_ = I*dx, J*dx
Xr, Yr = X0 + s*np.cos(th) - n_*np.sin(th), Y0 + s*np.sin(th) + n_*np.cos(th)
Z = np.zeros_like(s)
res = {}
for name, grid in [("cart", None), ("curv_id", (s, n_, Z)), ("rot30", (Xr, Yr, Z))]:
    ctrl = base + ("G 50 2\n" if grid is not None else "")
    d = os.path.join(S3, "runs/t4_"+name+var)
    rc = p1lib.run(d, control, ctrl, grid=grid, steps=steps)
    st = final_state(d)
    lines = p1lib.log_lines(d, r"^Inflow:|^Outflow:")
    qi = float(lines[-2].split()[1]) if len(lines)>=2 else float("nan")
    qo = float(lines[-1].split()[1]) if len(lines)>=2 else float("nan")
    res[name] = st
    print("%-8s rc=%d t=%.3f Qin %.4f Qout %.4f  max |U| %.3f  eta range %.4f %.4f" % (name, rc, st["simtime"], qi, qo, np.abs(st["U"]).max(), st["eta"].min(), st["eta"].max()))
for name, ref in [("curv_id","cart"), ("rot30","cart"), ("rot30","curv_id")]:
    r, c = res[name], res[ref]
    if name == "rot30":
        U = r["U"]*np.cos(th) + r["V"]*np.sin(th)
        V = -r["U"]*np.sin(th) + r["V"]*np.cos(th)
    else:
        U, V = r["U"], r["V"]
    print("%-8s vs %-7s: eta %.1e  U %.1e  V %.1e  W %.1e  P %.1e" % (name, ref, np.abs(r["eta"]-c["eta"]).max(), np.abs(U-c["U"]).max(), np.abs(V-c["V"]).max(), np.abs(r["W"]-c["W"]).max(), np.abs(r["P"]-c["P"]).max()))
