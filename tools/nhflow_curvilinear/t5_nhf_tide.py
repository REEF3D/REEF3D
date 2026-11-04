# P3 T5: oblique tide (direction 30 deg, period 600 s) through a 3 km basin with open edges on all four
# sides (Riemann at x- and y-, Flather at x+ and y+; regression case nhflow_3d_tide_oblique):
# Cartesian, curvilinear on the same nodes, turned by 30 deg at UTM coordinates (tide direction and
# origin turned with it). Per cell index, velocities turned back. XCTRL: extra ctrl lines (';')
import sys, numpy as np
from p3lib import *
import p1lib
steps = int(sys.argv[1]) if len(sys.argv) > 1 else 200
nx, ny, nz, dx, h = 60, 60, 3, 50.0, 10.0
L = nx*dx
control = "C 11 2\nC 12 2\nC 13 2\nC 14 2\nC 15 21\nC 16 3\n\nB 2 %d %d %d\nB 10 0.0 %g 0.0 %g 0.0 1.0\n\nM 10 1\nM 20 2\n" % (nx, ny, nz, L, L)
edges = "B 520 1 3 0\nB 521 1 1 0.0 0.0 0.0\nB 523 1 1\nB 520 3 3 0\nB 521 3 3 0.0 0.0 0.0\nB 523 3 1\nB 520 2 4 0\nB 521 2 2 0.0 0.0 0.0\nB 523 2 1\nB 520 4 4 0\nB 521 4 4 0.0 0.0 0.0\nB 523 4 1\n"
def ctrl(dirdeg, x0, y0):
    return (os.environ.get("XCTRL","").replace(";","\n") + "A 10 5\nA 520 %s\nB 90 1\nB 92 0\nB 510 1 1 %.1f 600.0\nB 511 1 0.25 600.0 90.0\nB 515 1 %.3f %.3f 0.0\n" % (os.environ.get("A520","2"), dirdeg, x0, y0)
            + edges + "F 60 %g\nN 47 0.5\nN 48 0\nN 49 2.0\nM 10 1\nP 10 0\nW 22 -9.81\n" % h)
X0, Y0 = 2.7e5, 7.04e6
th = np.radians(30.0)
I, J = np.meshgrid(np.arange(nx+1), np.arange(ny+1), indexing="ij")
s, n_ = I*dx, J*dx
Xr, Yr = X0 + s*np.cos(th) - n_*np.sin(th), Y0 + s*np.sin(th) + n_*np.cos(th)
Z = np.zeros_like(s)
res = {}
for name, grid, c in [("cart", None, ctrl(30.0, 0.0, 0.0)), ("curv_id", (s, n_, Z), ctrl(30.0, 0.0, 0.0)+"G 50 2\n"), ("rot30", (Xr, Yr, Z), ctrl(60.0, X0, Y0)+"G 50 2\n")]:
    d = os.path.join(S3, "runs/t5_"+name)
    rc = p1lib.run(d, control, c, grid=grid, steps=steps)
    st = final_state(d)
    res[name] = st
    print("%-8s rc=%d t=%.1f  eta range %.4f %.4f  max |U| %.4f |V| %.4f" % (name, rc, st["simtime"], st["eta"].min(), st["eta"].max(), np.abs(st["U"]).max(), np.abs(st["V"]).max()))
    for l in p1lib.log_lines(d, "!!!"):
        print("   ", l)
for name, ref in [("curv_id","cart"), ("rot30","cart"), ("rot30","curv_id")]:
    r, c = res[name], res[ref]
    if name == "rot30":
        U = r["U"]*np.cos(th) + r["V"]*np.sin(th); V = -r["U"]*np.sin(th) + r["V"]*np.cos(th)
    else:
        U, V = r["U"], r["V"]
    print("%-8s vs %-7s: eta %.1e  U %.1e  V %.1e  W %.1e  P %.1e" % (name, ref, np.abs(r["eta"]-c["eta"]).max(), np.abs(U-c["U"]).max(), np.abs(V-c["V"]).max(), np.abs(r["W"]-c["W"]).max(), np.abs(r["P"]-c["P"]).max()))
