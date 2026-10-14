# P5a T1: bed load over a flat erodible bed, flow from a collapsing water column (2D, radial): Cartesian,
# curvilinear on the same nodes, turned by 90 deg at (1000, 2000) m: bed and surface per cell index.
# environment: SOLVER nhflow|sflow, A520 (NHFLOW), XCTRL extra ctrl lines (;), NP ranks
import sys, numpy as np
from p5lib import *
steps = int(sys.argv[1]) if len(sys.argv) > 1 else 40
nx, ny, nz, dx = 32, 24, 6, 0.125
L, W = nx*dx, ny*dx
np_ = int(os.environ.get("NP", "1"))
solver = os.environ.get("SOLVER", "nhflow")
ctl = control(nx, ny, nz, L, W, 0.5, np_=np_)
if solver == "nhflow":
    base = "A 10 5\nA 520 %s\nN 47 0.3\n" % os.environ.get("A520", "0")
else:
    base = "A 10 2\nA 210 3\nA 211 4\nA 220 %s\nA 240 1\nN 47 0.3\n" % os.environ.get("A220", "0")
sed = "S 10 1\nS 11 2\nS 20 0.0005\nS 15 3\nS 17 1.0\nS 35 20.0\nS 41 1\nS 43 1\nS 44 1\n"
base += "N 48 0\nN 49 0.01\nF 60 0.2\nM 10 %d\nP 10 0\nW 22 -9.81\n" % np_ + sed + os.environ.get("XCTRL", "").replace(";", "\n")
box = (1.0, 1.75, 1.0, 2.0)
X0, Y0 = 1000.0, 2000.0
res = {}
for name, rot in [("cart", None), ("curv_id", "id"), ("rot90", "rot90")]:
    if rot is None:
        grid = None
        b = box
    else:
        X, Y = nodes(nx, ny, dx, rot, X0, Y0)
        grid = (X, Y, np.zeros_like(X))
        b = box_phys(box, rot, X0, Y0)
    ctrl = base + "F 72 %g %g %g %g 0.35\n" % b + ("G 50 2\n" if grid is not None else "")
    d = os.path.join(S5, "runs/t1_%s_%s" % (solver, name))
    rc = run(d, ctl, ctrl, grid=grid, np_=np_, steps=steps)
    st = final_state(d)
    res[name] = st
    db = st["bed"]
    print("%-8s rc=%d t=%.3f  bed min %.3e max %.3e  sum %.3e  eta max %.4f" % (name, rc, st["simtime"], db.min(), db.max(), db.sum(), st["eta"].max()),
          log_lines(d, "!!!")[:2])
for name, ref in [("curv_id","cart"), ("rot90","curv_id"), ("rot90","cart")]:
    r, c = res[name], res[ref]
    scale = max(np.abs(c["bed"]).max(), 1e-30)
    print("%-8s vs %-7s: bed %.1e (rel %.1e)  eta %.1e  WL %.1e" % (name, ref, np.abs(r["bed"]-c["bed"]).max(), np.abs(r["bed"]-c["bed"]).max()/scale,
          np.abs(r["eta"]-c["eta"]).max(), np.abs(r["WL"]-c["WL"]).max()))
