# P5c T1: 6DOF bodies on a curvilinear grid (G 50 2) against the Cartesian grid: the Cartesian run, the curvilinear
# grid on the same nodes (curv_id) and the grid turned by 90 deg about (X0, Y0) (rot90). The body is the same
# box (X 110) about its centre of gravity; in the turned run its centre is turned and its initial yaw is 90 deg
# more (X 101: the body frame stays that of the box), the initial linear velocity (X 102, inertial) and the
# prescribed velocity (X 210) are turned, the angular velocity (X 103) is body-fixed; the prescribed degrees of
# freedom (X 11) are inertial: surge and sway both prescribed. The box pressure of the ship (X 400 2) is aligned with
# x (not turned with the grid): the turned runs compare with CASE shipwave_stl (X 400 10, the pressure of the hull).
#   curv_id vs cart: per cell index, eta, U, V, W (NHFLOW) / eta, WL (SFLOW): round-off for NHFLOW; SFLOW builds its
#                    2D body level set differently on the curvilinear grid (exact distance to the waterline instead
#                    of rays along the grid lines and re-initialisation): close, not round-off
#   rot90 vs curv_id: per cell index (velocities turned back): round-off
#   motions: the centre of gravity, the angles and the velocities of the body (6 digits of the output files)
# fixed time step (N 48 0, N 49): NHFLOW takes the adaptive time step from the speed on curvilinear grids
# rot30: the WENO weights (A 514 4, SFLOW A 211) act on the Cartesian components and the default direct forcing
# X 15 1 weights each Cartesian component by |n_x|, |n_y|: not invariant under a turn other than 90 deg; rot30 with
# XCTRL "A 514 0;A 524 1;X 15 2" (NHFLOW) or "A 211 0" (SFLOW)
# environment: CASE floating|towed|thrust|propeller|shipwave|shipwave_stl (NHFLOW) / towed|shipwave|shipwave_stl (SFLOW), SOLVER,
# NP, XCTRL (;), NORUN
import sys, os, glob, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../sediment_curvilinear"))
from p5lib import control, nodes, run5, final_state, log_lines

T6 = os.path.dirname(os.path.abspath(__file__))
steps = int(sys.argv[1]) if len(sys.argv) > 1 else 40
solver = os.environ.get("SOLVER", "nhflow")
case = os.environ.get("CASE", "floating")
np_ = int(os.environ.get("NP", "1"))
X0, Y0 = 1000.0, 2000.0

nx, ny, nz, dx = 80, 40, 10, 0.1
L, W, h = nx*dx, ny*dx, 0.5
xc, yc = 3.0, 2.0                      # centre of the box (x, y) in the unturned frame
bx, by, zb0, zb1 = 0.8, 0.4, 0.4, 0.6   # box size and z extent
yaw0 = float(os.environ.get("YAW", "10.0"))
lin = tuple(float(v) for v in os.environ.get("LIN", "0.10,0.04,0.0").split(","))   # X 102 (inertial)
ang = tuple(float(v) for v in os.environ.get("ANG", "0.15,0.10,0.05").split(","))  # X 103 (body-fixed)
files = {}
bc = "C 11 21\nC 12 21\nC 13 21\nC 14 21\nC 15 21\nC 16 3\n"

if solver == "nhflow":
    flow = "A 10 5\nN 47 0.3\nN 48 0\nN 49 0.01\nN 41 60.0\nM 10 %d\nF 60 %g\nW 22 -9.81\n" % (np_, h)
    if case == "floating":
        body = "X 10 1\nX 11 1 1 1 1 1 1\nX 21 500.0\n"
    elif case == "towed":
        body = "X 10 2\nX 11 2 2 0 0 0 0\n"
    elif case in ("thrust", "propeller"):
        body = "X 10 1\nX 11 1 1 1 1 1 1\nX 21 500.0\nX 350 1\n"
        bx, by, zb0, zb1 = 1.0, 0.3, 0.35, 0.65
        files["ship.dat"] = ("form_factor 0.1\nthrust 1.0 -0.5 -0.03\n" if case == "thrust" else
                             "form_factor 0.1\npropeller -0.6 0.0 -0.07 0.12\npropeller_hub 0.2\npropeller_kt 0.33 -0.25 -0.08\n"
                             "propeller_kq 0.045 -0.03 -0.01\npropeller_rps 6.0\npropeller_inflow wake 0.2\n")
        lin, ang = (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)
    elif case == "shipwave":
        body = "X 10 3\nX 11 2 2 0 0 0 0\nX 400 2\nX 401 200.0 2.0 16.0 16.0\n"
        bx, by, zb0, zb1 = 1.0, 0.5, 0.45, 0.55   # through the free surface (F 60 0.5)
        yaw0 = 0.0
    elif case == "shipwave_stl":
        body = "X 10 3\nX 11 2 2 0 0 0 0\nX 400 10\nX 401 200.0 2.0 16.0 16.0\n"
        bx, by, zb0, zb1 = 1.0, 0.5, 0.45, 0.55
else:
    flow = "A 10 2\nA 210 3\nA 211 4\nA 220 0\nA 240 1\nN 47 0.3\nN 48 0\nN 49 0.01\nN 41 60.0\nM 10 %d\nF 60 0.4\nW 22 -9.81\n" % np_
    h = 0.4
    bx, by, zb0, zb1 = 1.5, 0.5, 0.3, 0.5
    if case == "towed":
        body = "X 10 2\nX 11 2 2 0 0 0 0\n"
    elif case == "shipwave":
        body = "X 10 3\nX 11 2 2 0 0 0 0\nX 400 2\nX 401 300.0 2.0 16.0 16.0\n"
        yaw0 = 0.0
    elif case == "shipwave_stl":
        body = "X 10 3\nX 11 2 2 0 0 0 0\nX 400 10\n"

# the turned runs: (angle, X0, Y0); rot90 about (1000, 2000) (exact), rot30 at UTM coordinates (ROT30=1) or about (1000, 2000) (ROT30=small)
TURN = {"rot90": (90.0, X0, Y0), "rot30": (30.0, 2.7e5, 7.04e6) if os.environ.get("ROT30") != "small" else (30.0, X0, Y0)}

def turn(v, rot):
    if rot not in TURN:
        return v
    if rot == "rot90":
        return (-v[1], v[0], v[2])
    t = np.radians(TURN[rot][0])
    return (v[0]*np.cos(t) - v[1]*np.sin(t), v[0]*np.sin(t) + v[1]*np.cos(t), v[2])

def place(x, y, rot):
    if rot not in TURN:
        return x, y
    th, ox, oy = TURN[rot]
    if rot == "rot90":
        return ox - y, oy + x
    t = np.radians(th)
    return ox + x*np.cos(t) - y*np.sin(t), oy + x*np.sin(t) + y*np.cos(t)

def body_ctrl(rot):
    cx, cy = place(xc, yc, rot)
    c = "X 110 %.17g %.17g %.17g %.17g %g %g\n" % (cx-0.5*bx, cx+0.5*bx, cy-0.5*by, cy+0.5*by, zb0, zb1)
    yaw = yaw0 + (TURN[rot][0] if rot in TURN else 0.0)
    # FIXCG=1: centre of gravity and inertia given (X 23, X 24): the comparison without the round-off of the
    # surface integrals (about the origin on the Cartesian grid, about the body on the curvilinear grid)
    if os.environ.get("FIXCG"):
        c += "X 23 %.17g %.17g %g\nX 24 2.0 4.0 5.0\n" % (cx, cy, 0.5*(zb0+zb1))
    if yaw != 0.0:
        c += "X 101 0.0 0.0 %.17g\n" % yaw
    if case in ("floating",):
        c += "X 102 %.17g %.17g %.17g\nX 103 %g %g %g\n" % (turn(lin, rot) + ang)
    if case in ("towed", "shipwave", "shipwave_stl"):
        c += "X 210 %.17g %.17g %.17g\n" % turn((0.3, 0.05 if case == "towed" else 0.0, 0.0), rot)
    return c

ctl = control(nx, ny, nz, L, W, 1.0, bc=bc, np_=np_)
res, mot = {}, {}
RUNS = [("cart", None), ("curv_id", "id"), ("rot90", "rot90")] + ([("rot30", "rot30")] if os.environ.get("ROT30") else [])
for name, rot in RUNS:
    grid = None
    ctrl = flow + body + body_ctrl(rot) + os.environ.get("XCTRL", "").replace(";", "\n") + "\n"
    if rot is not None:
        X, Y = nodes(nx, ny, dx, "rot", TURN[rot][1], TURN[rot][2], TURN[rot][0]) if rot == "rot30" else nodes(nx, ny, dx, rot, X0, Y0)
        grid = (X, Y, np.zeros_like(X))
        ctrl += "G 50 2\n"
    d = os.path.join(T6, "runs/t1_%s_%s_%s_np%d" % (solver, case, name, np_))
    rc = 0
    if not os.environ.get("NORUN"):
        rc = run5(d, ctl, ctrl, grid=grid, np_=np_, steps=steps, files=files)
    st = [final_state(d, r) for r in range(np_)]
    res[name] = st
    m = {}
    for kind in ("position", "velocity", "forces"):
        f = glob.glob(os.path.join(d, "REEF3D_*_6DOF", "REEF3D_6DOF_%s_0.dat" % kind))
        if f:
            a = np.loadtxt(f[0], skiprows=1, ndmin=2)
            m[kind] = a[-1] if len(a) else None
    mot[name] = m
    print("%-8s rc=%d t=%.3f eta max %.4e" % (name, rc, st[0]["simtime"], max(np.abs(s["eta"]).max() for s in st)), log_lines(d, "!!!")[:2])

fields = ["eta", "U", "V", "W"] if solver == "nhflow" else ["eta", "WL"]
def cmp(a, b, back=False):
    out = []
    for f in fields:
        e = 0.0
        for r in range(np_):
            fa, fb = res[a][r][f], res[b][r][f]
            if back and f == "U":
                fa = res[a][r]["V"]
            if back and f == "V":
                fa = -res[a][r]["U"]
            e = max(e, np.abs(fa-fb).max())
        out.append("%s %.1e" % (f, e))
    return "  ".join(out)
print("curv_id vs cart   :", cmp("curv_id", "cart"))
print("rot90   vs curv_id:", cmp("rot90", "curv_id", back=True))
if "rot30" in res:
    t = np.radians(30.0)
    for r in range(np_):
        U, V = res["rot30"][r]["U" if solver=="nhflow" else "eta"], None
    if solver == "nhflow":
        e = {}
        for r in range(np_):
            a, b = res["rot30"][r], res["curv_id"][r]
            Ub = a["U"]*np.cos(t) + a["V"]*np.sin(t)
            Vb = -a["U"]*np.sin(t) + a["V"]*np.cos(t)
            e["eta"] = max(e.get("eta",0), np.abs(a["eta"]-b["eta"]).max())
            e["U"] = max(e.get("U",0), np.abs(Ub-b["U"]).max())
            e["V"] = max(e.get("V",0), np.abs(Vb-b["V"]).max())
        print("rot30   vs curv_id:", "  ".join("%s %.1e" % kv for kv in e.items()))
    else:
        print("rot30   vs curv_id: eta %.1e" % max(np.abs(res["rot30"][r]["eta"]-res["curv_id"][r]["eta"]).max() for r in range(np_)))
for name, _ in RUNS:
    m = mot[name]
    if "position" in m and m["position"] is not None:
        pos = m["position"].copy()
        if name in TURN:
            th, ox, oy = TURN[name]
            t = np.radians(th)
            x, y = pos[1]-ox, pos[2]-oy
            pos[1], pos[2] = x*np.cos(t) + y*np.sin(t), -x*np.sin(t) + y*np.cos(t)
            pos[6] -= th
        if name in TURN:   # 6 significant digits: the turned centre (1e3, 1e6) is not resolved, only zg and the angles
            print("%-8s position t=%.3f zg %.6f phi %.4f theta %.4f psi %.4f" % ((name, pos[0]) + tuple(pos[3:7])))
        else:
            print("%-8s position t=%.3f xg %.6f yg %.6f zg %.6f phi %.4f theta %.4f psi %.4f" % ((name,) + tuple(pos[:7])))
