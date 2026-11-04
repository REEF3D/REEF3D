# P3 T6: several ranks. (a) the dam break of T1 (A 520 2) on 2 and 4 ranks: curvilinear on the
# Cartesian nodes against the Cartesian run per rank to round-off; (b) hydrostatic (A 520 0) on the
# grid turned by 30 deg: 2 and 4 ranks against one rank, bit for bit (sorted values; the
# non-hydrostatic runs differ between rank counts by the partition dependence of the pressure solve,
# as on Cartesian grids)
import sys, numpy as np
from p3lib import *
import p1lib
steps = int(sys.argv[1]) if len(sys.argv) > 1 else 20
n, nz, dx = 40, 8, 0.05
I, J = np.meshgrid(np.arange(n+1), np.arange(n+1), indexing="ij")
X0, Y0 = I*dx, J*dx; Z = np.zeros_like(X0)
def control(np_):
    return "C 11 21\nC 12 21\nC 13 21\nC 14 21\nC 15 21\nC 16 3\n\nB 2 %d %d %d\nB 10 0.0 %g 0.0 %g 0.0 0.5\n\nM 10 %d\n" % (n, n, nz, n*dx, n*dx, np_)
for np_ in [2, 4]:
    ctrl = "A 10 5\nA 520 2\nF 60 0.5\nF 72 0.91 1.11 0.91 1.11 0.54\nN 47 0.5\nM 10 %d\nP 10 0\nW 22 -9.81\n" % np_
    res = {}
    for name, g in [("cart", None), ("curv", (X0, Y0, Z))]:
        d = os.path.join(S3, "runs/t6_%s_%d" % (name, np_))
        rc = p1lib.run(d, control(np_), ctrl + ("G 50 2\n" if g else ""), grid=g, np_=np_, steps=steps)
        res[name] = [final_state(d, r) for r in range(np_)]
    print("(a) np %d curv vs cart: " % np_ + "  ".join("%s %.1e" % (f, max(np.abs(res["curv"][r][f]-res["cart"][r][f]).max() for r in range(np_))) for f in ["eta","U","V","W","P"]))
th = np.radians(30.0)
Xr = 1000.0 + X0*np.cos(th) - Y0*np.sin(th); Yr = 2000.0 + X0*np.sin(th) + Y0*np.cos(th)
cx, cy = 1000.0 + np.cos(th) - np.sin(th), 2000.0 + np.sin(th) + np.cos(th)
ref = None
for np_ in [1, 2, 4]:
    ctrl = "A 10 5\nA 520 0\nF 60 0.5\nF 72 %.4f %.4f %.4f %.4f 0.54\nN 47 0.5\nN 48 0\nN 49 0.005\nM 10 %d\nP 10 0\nW 22 -9.81\nG 50 2\n" % (cx-0.12, cx+0.12, cy-0.12, cy+0.12, np_)
    d = os.path.join(S3, "runs/t6_rot30_np%d" % np_)
    rc = p1lib.run(d, control(np_), ctrl, grid=(Xr, Yr, Z), np_=np_, steps=steps)
    v = {f: np.sort(np.concatenate([final_state(d, r)[f] for r in range(np_)])) for f in ["eta","U","V","W"]}
    if ref is None:
        ref = v
        continue
    print("(b) np %d vs np 1 (A 520 0, turned by 30 deg): " % np_ + "  ".join("%s %.1e" % (f, np.abs(v[f]-ref[f]).max()) for f in v))
