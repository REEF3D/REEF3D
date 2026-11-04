# P5c T3: 6DOF bodies on a multi-block grid (G 50 1, DIVEMesh R 20 1): the skewed wavy grid of P5a T5 as one block
# and as the four quadrants with their index frames turned by 0, 90, 180, 270 deg; a box over the vertex of the four
# blocks (the hull lies in all of them). Per global cell (level, I, J): eta, U, V (NHFLOW: the bottom layer); the
# motion of the body.
#   four vs one: the same solution to round-off (the loads of every hull triangle from the block that holds its
#                centroid, the hull cast on the grid of every block; the surface remeshing of the hull, X 185, and
#                the bucket order of the triangles do not depend on the blocks)
# environment: SOLVER nhflow|sflow, CASE floating|towed (NHFLOW) / towed|shipwave_stl (SFLOW), NP (ranks of the
# four-block run), LEV (G 1 with the zone around the hull G 12), NOTURN (the blocks not turned), SHIFTX, SHIFTY (the
# box moved from the vertex), ROTONE (one block with its index frame turned by 90 deg instead of the four), TOP (the cells of the largest differences of eta), XCTRL (;), NORUN
import sys, os, glob, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../sediment_curvilinear"))
from p5lib import amr_grids2, read_state, log_lines
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../nhflow_multiblock"))
from p4lib import run4

T6 = os.path.dirname(os.path.abspath(__file__))
steps = int(sys.argv[1]) if len(sys.argv) > 1 else 30
solver = os.environ.get("SOLVER", "nhflow")
case = os.environ.get("CASE", "floating" if solver == "nhflow" else "towed")
np_ = int(os.environ.get("NP", "1"))
LEV = int(os.environ.get("LEV", "0"))

n, nz, dx = 40, 6, 0.1
L = n*dx
h = n//2
I, J = np.meshgrid(np.arange(n+1), np.arange(n+1), indexing="ij")
x, y = I*dx, J*dx
X = x + 0.15*y + 0.04*np.sin(2*np.pi*x/L)*np.sin(np.pi*y/L)
Y = y + 0.05*x + 0.03*np.sin(np.pi*x/L)*np.sin(2*np.pi*y/L)
Z = np.zeros_like(X)
xv, yv = X[h, h], Y[h, h]           # the vertex of the four blocks

def quad(i0, j0, k):
    return tuple(np.rot90(A[i0:i0+h+1, j0:j0+h+1], k=-k).copy() for A in (X, Y, Z))

def rot(m, k, a, b):
    q = np.rot90(np.arange(m*m).reshape(m, m), k=-k)[a, b]
    return q//m, q % m

offs = [(0, 0, 0), (h, 0, 1), (0, h, 2), (h, h, 3)] if not os.environ.get("NOTURN") else [(0, 0, 0), (h, 0, 0), (0, h, 0), (h, h, 0)]
blocks = [quad(*o) for o in offs]
control = "B 1 %g\nB 10 0.0 %g 0.0 %g 0.0 %g\n\nM 10 NP\n\nR 20 1\n" % (dx, L, L, 1.0)

bx, by, zb0, zb1 = 0.8, 0.4, 0.4, 0.6
cx, cy = xv + 0.05 + float(os.environ.get("SHIFTX", "0")), yv + 0.03 + float(os.environ.get("SHIFTY", "0"))
if solver == "nhflow":
    hw = 0.5
    flow = "A 10 5\n"
    if case == "floating":
        body = "X 10 1\nX 11 1 1 1 1 1 1\nX 21 500.0\nX 102 0.10 0.04 0.0\nX 103 0.15 0.10 0.05\n"
    else:
        body = "X 10 2\nX 11 2 2 0 0 0 0\nX 210 0.3 0.05 0.0\n"
else:
    hw = 0.4
    flow = "A 10 2\nA 209 0\nA 210 3\nA 211 4\nA 220 0\nA 240 1\n"
    bx, by, zb0, zb1 = 1.2, 0.5, 0.3, 0.5
    if case == "towed":
        body = "X 10 2\nX 11 2 2 0 0 0 0\nX 210 0.3 0.05 0.0\n"
    else:
        body = "X 10 3\nX 11 2 2 0 0 0 0\nX 400 10\nX 210 0.3 0.0 0.0\n"
body += "X 110 %.17g %.17g %.17g %.17g %g %g\nX 101 0.0 0.0 10.0\n" % (cx-0.5*bx, cx+0.5*bx, cy-0.5*by, cy+0.5*by, zb0, zb1)
amr = "G 1 %d\nG 2 4\nG 3 2\nG 4 8\nG 12 0.3\n" % LEV if LEV > 0 else ""
ctrl = (flow + "F 60 %g\nG 50 1\nG 52 1\nN 47 0.3\nN 48 0\nN 49 0.01\nN 41 60.0\nM 10 NP\nP 10 0\nW 22 -9.81\n" % hw
        + body + amr + os.environ.get("XCTRL", "").replace(";", "\n") + "\n")

def grids(d, npr):
    if os.path.exists(os.path.join(d, "amr", "amr_r0.bin")):
        return amr_grids2(d, npr)
    out = []
    for r in range(npr):
        st = read_state(sorted(glob.glob(os.path.join(d, "dump", "state_*_r%d.bin" % r)))[-1])
        G = {"rank": r, "id": 0, "lev": 0, "blk": 0, "I0": 0, "J0": 0, "nx": n, "ny": n}
        G["eta"] = np.asarray(st["eta"]).reshape(n, n)
        G["U"] = G["V"] = np.full((n, n), np.nan)     # SFLOW: face fluxes only in the regression dump
        out.append(G)
    return out

def cell_key(G4, g, a, b, kind):
    lv = g["lev"]
    m = h << lv
    if kind == "one":
        return (lv, g["I0"]+a, g["J0"]+b)
    if kind == "rotone":
        p, q = rot(n << lv, 1, g["I0"]+a, g["J0"]+b)
        return (lv, p, q)
    O = {gg["blk"]: (gg["I0"], gg["J0"]) for gg in G4 if gg["lev"] == 0 and gg["blk"] > 0}
    O[0] = (0, 0)
    A0 = g["I0"]+a - (O[g["blk"]][0] << lv)
    B0 = g["J0"]+b - (O[g["blk"]][1] << lv)
    if A0 < 0 or B0 < 0 or A0 >= m or B0 >= m:
        return None
    i0, j0, k = offs[g["blk"]]
    p, q = rot(m, k, A0, B0)
    return (lv, (i0 << lv) + p, (j0 << lv) + q)

T, mot = {}, {}
RUNS = [("one", [(X, Y, Z)], 1, "one"), ("four", blocks, np_, "four")]
if os.environ.get("ROTONE"):
    RUNS[1] = ("four", [tuple(np.rot90(A, k=-1).copy() for A in (X, Y, Z))], np_, "rotone")
for name, bl, npr, kind in RUNS:
    d = os.path.join(T6, "runs/t3_%s_%s_%s_L%d_np%d" % (solver, case, name, LEV, npr))
    c = ctrl + ("G 1 1\n" if len(bl) == 1 and LEV == 0 else "")
    rc = 0
    if not os.environ.get("NORUN"):
        rc = run4(d, control.replace("NP", str(npr)), c.replace("NP", str(npr)), bl, np_=npr, steps=steps)
    G = grids(d, npr)
    t = {}
    for g in G:
        for a in range(g["nx"]):
            for b in range(g["ny"]):
                key = cell_key(G, g, a, b, kind)
                if key is not None:
                    t[key] = (g["eta"][a, b], g["U"][a, b], g["V"][a, b])
    T[name] = t
    f = glob.glob(os.path.join(d, "REEF3D_*_6DOF", "REEF3D_6DOF_position_0.dat"))
    a = np.loadtxt(f[0], skiprows=1, ndmin=2) if f else np.zeros((0, 7))
    mot[name] = a[-1] if len(a) else None
    levs = sorted({k[0] for k in t})
    print("%-4s rc=%d grids %d levels %s cells %s" % (name, rc, len(G), levs, [sum(1 for k in t if k[0] == l) for l in levs]),
          log_lines(d, "!!!")[:3])

for lv in sorted({k[0] for k in T["one"]}):
    ka = {k for k in T["one"] if k[0] == lv}
    kb = {k for k in T["four"] if k[0] == lv}
    com = ka & kb
    e = [max([abs(T["four"][k][m]-T["one"][k][m]) for k in com] or [float("nan")]) for m in range(3)]
    print("four vs one (level %d, %d cells%s): eta %.1e U %.1e V %.1e" % (lv, len(com), ", %d/%d not common" % (len(ka-kb), len(kb-ka)) if ka != kb else "", *e))
for name in ("one", "four"):
    m = mot[name]
    if m is not None:
        print("%-4s position t=%.3f xg %.6f yg %.6f zg %.6f phi %.4f theta %.4f psi %.4f" % ((name,) + tuple(m[:7])))
if os.environ.get("TOP"):
    e = sorted(((abs(T["four"][k][0]-T["one"][k][0]), k) for k in T["one"] if k in T["four"]), reverse=True)[:12]
    print("largest |deta| at", ["%s %.1e" % (k, v) for v, k in e])
