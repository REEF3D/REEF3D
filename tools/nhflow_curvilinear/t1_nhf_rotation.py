# P3 T1: NHFLOW dam break (F 72 hump) in a closed basin: Cartesian, curvilinear on the same nodes,
# curvilinear on the grid turned by 90 deg and shifted: per cell index the same to round-off
import sys, numpy as np
from p3lib import *
import p1lib
steps = int(sys.argv[1]) if len(sys.argv) > 1 else 40
n, nz, dx = 40, 8, 0.05
control = "C 11 21\nC 12 21\nC 13 21\nC 14 21\nC 15 21\nC 16 3\n\nB 2 %d %d %d\nB 10 0.0 %g 0.0 %g 0.0 0.5\n\nM 10 1\n" % (n, n, nz, n*dx, n*dx)
ctrl = os.environ.get("CTRL", "A 10 5\nF 60 0.5\nF 72 0.91 1.11 0.91 1.11 0.54\nN 47 0.5\nM 10 1\nP 10 0\nW 22 -9.81\n")
I, J = np.meshgrid(np.arange(n+1), np.arange(n+1), indexing="ij")
X0, Y0 = I*dx, J*dx
Z = np.zeros_like(X0)
Xr = float(os.environ.get("OX","1000.0")) - (Y0-1.0)
Yr = float(os.environ.get("OY","2000.0")) + (X0-1.0)
# turned by 180 deg: the walls stay walls of the same index direction (the Cartesian NHFLOW walls are
# not isotropic: x-walls set the ghost velocities to 0, y-walls are slip walls), 90 deg turns them
Xh = float(os.environ.get("OX","1000.0")) - (X0-1.0)
Yh = float(os.environ.get("OY","2000.0")) - (Y0-1.0)
res = {}
for name, g, grid in [("cart", 0, None), ("curv_id", 2, (X0, Y0, Z)), ("curv_180", 2, (Xh, Yh, Z)), ("curv_rot", 2, (Xr, Yr, Z))]:
    c = ctrl + ("G 50 %d\n" % g if g else "")
    if name in ("curv_rot", "curv_180"):
        ox, oy = float(os.environ.get("OX","1000.0")), float(os.environ.get("OY","2000.0"))
        c = c.replace("F 72 0.91 1.11 0.91 1.11", "F 72 %.6f %.6f %.6f %.6f" % (ox-0.09, ox+0.11, oy-0.09, oy+0.11))
    d = os.path.join(S3, "runs/t1_"+name)
    rc = p1lib.run(d, control, c, grid=grid, steps=steps)
    st = final_state(d)
    res[name] = st
    print("%-9s rc=%d count=%d t=%.6f fields %s" % (name, rc, st["count"], st["simtime"], sorted(k for k in st if k not in ("count","simtime"))))
for name in ["curv_id", "curv_180", "curv_rot"]:
    r, c = res[name], res["cart"]
    out = []
    for f in ["eta", "WL", "U", "V", "W", "P"]:
        a, b = r[f], c[f]
        if name == "curv_rot" and f == "U": b = -c["V"]
        if name == "curv_rot" and f == "V": b = c["U"]
        if name == "curv_180" and f in ("U", "V"): b = -c[f]
        out.append("%s %.1e" % (f, np.abs(a-b).max()))
    print("%-9s vs cart (max |eta| %.3f, |U| %.3f): %s" % (name, np.abs(c["eta"]).max(), np.abs(c["U"]).max(), "  ".join(out)))
