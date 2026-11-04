# P5c T2: 6DOF bodies with mesh refinement (G 1) on a curvilinear grid (G 50 2): the box of T1 with a refinement
# zone around the hull (G 12, the patches follow the body, regrid every G 2 steps), Cartesian, curvilinear on the
# same nodes (curv_id) and turned by 90 deg about (1000, 2000) (rot90). Per global cell (level, I, J) of all grids of
# all ranks (the AMR dump): eta, U, V (NHFLOW: the bottom layer); the motion of the body.
#   curv_id vs cart:  NHFLOW round-off (the surface remeshing of the hull, X 185, breaks ties differently: 1e-5);
#                     SFLOW the level set of the waterline is built differently (as T1)
#   rot90 vs curv_id: round-off, the same patches
# environment: CASE floating|towed (NHFLOW) / towed|shipwave_stl (SFLOW), SOLVER, NP, LEV (G 1), WAKE (G 13 length:
# the zone along the heading of the body), PLACE (G 40 1: NHFLOW patches placed for the load), XCTRL (;), NORUN
import sys, os, glob, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../sediment_curvilinear"))
from p5lib import control, nodes, run5, amr_grids2, log_lines

T6 = os.path.dirname(os.path.abspath(__file__))
steps = int(sys.argv[1]) if len(sys.argv) > 1 else 24
solver = os.environ.get("SOLVER", "nhflow")
case = os.environ.get("CASE", "floating" if solver == "nhflow" else "towed")
np_ = int(os.environ.get("NP", "1"))
LEV = int(os.environ.get("LEV", "1"))
X0, Y0 = 1000.0, 2000.0

nx, ny, nz, dx = 80, 40, 6, 0.1
xc, yc = 3.0, 2.0
bx, by, zb0, zb1 = 0.8, 0.4, 0.4, 0.6
yaw0 = float(os.environ.get("YAW", "10.0"))
bc = "C 11 21\nC 12 21\nC 13 21\nC 14 21\nC 15 21\nC 16 3\n"
amr = "G 1 %d\nG 2 4\nG 3 2\nG 4 8\nG 12 0.3\n" % LEV
if os.environ.get("WAKE"):
    amr += "G 13 %s 19.47\n" % os.environ["WAKE"]
if os.environ.get("PLACE"):
    amr += "G 40 1\n"

if solver == "nhflow":
    h = 0.5
    flow = "A 10 5\nN 47 0.3\nN 48 0\nN 49 0.01\nN 41 60.0\nM 10 %d\nF 60 %g\nW 22 -9.81\n" % (np_, h)
    if case == "floating":
        body = "X 10 1\nX 11 1 1 1 1 1 1\nX 21 500.0\n"
    elif case == "towed":
        body = "X 10 2\nX 11 2 2 0 0 0 0\n"
else:
    h = 0.4
    flow = "A 10 2\nA 210 3\nA 211 4\nA 220 0\nA 240 1\nN 47 0.3\nN 48 0\nN 49 0.01\nN 41 60.0\nM 10 %d\nF 60 %g\nW 22 -9.81\n" % (np_, h)
    bx, by, zb0, zb1 = 1.5, 0.5, 0.3, 0.5
    if case == "towed":
        body = "X 10 2\nX 11 2 2 0 0 0 0\n"
    elif case == "shipwave_stl":
        body = "X 10 3\nX 11 2 2 0 0 0 0\nX 400 10\n"

def place(x, y, rot):
    return (X0 - y, Y0 + x) if rot == "rot90" else (x, y)

def turn(v, rot):
    return (-v[1], v[0], v[2]) if rot == "rot90" else v

def body_ctrl(rot):
    cx, cy = place(xc, yc, rot)
    c = "X 110 %.17g %.17g %.17g %.17g %g %g\n" % (cx-0.5*bx, cx+0.5*bx, cy-0.5*by, cy+0.5*by, zb0, zb1)
    yaw = yaw0 + (90.0 if rot == "rot90" else 0.0)
    if yaw != 0.0:
        c += "X 101 0.0 0.0 %.17g\n" % yaw
    if case == "floating":
        c += "X 102 %.17g %.17g %.17g\nX 103 0.15 0.10 0.05\n" % turn((0.10, 0.04, 0.0), rot)
    else:
        c += "X 210 %.17g %.17g %.17g\n" % turn((0.3, 0.05 if case == "towed" else 0.0, 0.0), rot)
    return c

ctl = control(nx, ny, nz, nx*dx, ny*dx, 1.0, bc=bc, np_=np_)
T, mot = {}, {}
RUNS = [("cart", None), ("curv_id", "id"), ("rot90", "rot90")]
for name, rot in RUNS:
    grid = None
    ctrl = flow + amr + body + body_ctrl(rot) + os.environ.get("XCTRL", "").replace(";", "\n") + "\n"
    if rot is not None:
        X, Y = nodes(nx, ny, dx, rot, X0, Y0)
        grid = (X, Y, np.zeros_like(X))
        ctrl += "G 50 2\n"
    d = os.path.join(T6, "runs/t2_%s_%s_%s_L%d_np%d" % (solver, case, name, LEV, np_))
    rc = 0
    if not os.environ.get("NORUN"):
        rc = run5(d, ctl, ctrl, grid=grid, np_=np_, steps=steps)
    G = amr_grids2(d, np_) if glob.glob(os.path.join(d, "amr", "amr_r*.bin")) else []
    t = {}
    for g in G:
        for a in range(g["nx"]):
            for b in range(g["ny"]):
                U, V = g["U"][a, b], g["V"][a, b]
                if rot == "rot90":
                    U, V = V, -U
                t[(g["lev"], g["I0"]+a, g["J0"]+b)] = (g["eta"][a, b], U, V)
    T[name] = t
    f = glob.glob(os.path.join(d, "REEF3D_*_6DOF", "REEF3D_6DOF_position_0.dat"))
    a = np.loadtxt(f[0], skiprows=1, ndmin=2) if f else np.zeros((0, 7))
    mot[name] = a[-1] if len(a) else None
    levs = sorted({k[0] for k in t})
    print("%-8s rc=%d grids %d levels %s cells %s" % (name, rc, len(G), levs, [sum(1 for k in t if k[0] == l) for l in levs]),
          log_lines(d, "!!!")[:2])

def cmp(a, b):
    out = []
    for lv in sorted({k[0] for k in T[b]}):
        ka = {k for k in T[a] if k[0] == lv}
        kb = {k for k in T[b] if k[0] == lv}
        com = ka & kb
        e = [max([abs(T[a][k][m]-T[b][k][m]) for k in com] or [float("nan")]) for m in range(3)]
        out.append("L%d (%d cells%s) eta %.1e U %.1e V %.1e" % (lv, len(com), ", %d/%d not common" % (len(ka-kb), len(kb-ka)) if ka != kb else "", *e))
    return "; ".join(out)
print("curv_id vs cart   :", cmp("curv_id", "cart"))
print("rot90   vs curv_id:", cmp("rot90", "curv_id"))
for name, _ in RUNS:
    m = mot[name]
    if m is not None:
        if name == "rot90":
            print("%-8s position t=%.3f zg %.6f phi %.4f theta %.4f psi %.4f" % (name, m[0], m[3], m[4], m[5], m[6]-90.0))
        else:
            print("%-8s position t=%.3f xg %.6f yg %.6f zg %.6f phi %.4f theta %.4f psi %.4f" % ((name,) + tuple(m[:7])))
