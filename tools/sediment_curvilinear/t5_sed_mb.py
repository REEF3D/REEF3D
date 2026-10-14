# P5a T5: sediment on multi-block grids (G 50 1, DIVEMesh R 20 1): the four quadrants of P4 T2 (index frames
# turned by 0, 90, 180, 270 deg, wavy skewed grid, bumpy bed), a dam break with bed load over the vertex:
# one block vs. four, per cell: bed, eta, U, V at level 0 (LEV > 0: and every fine cell); the bed volume change
# of both (sum of cell area times bed change; the initial bed from a run of 1 step with the morphological
# factor 0). The filter (S 100, default on) and the sand slide (XCTRL "S 90 1".."S 90 5", steep bump BUMP=0.2)
# act across the block interfaces.
# environment: SOLVER nhflow|sflow, A520, NP (ranks of the four-block run), LEV (G 1, static box over the vertex),
# BUMP (bump height, 0.06), XCTRL extra ctrl lines separated by ;, NOTURN=1 (the blocks not turned), ROTONE=1 (one
# block with its index frame turned by 90 deg instead of the four blocks), TOP=1 (the cells of the largest bed
# differences), NORUN=1 (evaluate the runs that are there)
import sys, glob, numpy as np
from p5lib import *
sys.path.insert(0, os.path.join(S5, "../nhflow_multiblock"))
from p4lib import run4

steps = int(sys.argv[1]) if len(sys.argv) > 1 else 40
np_ = int(os.environ.get("NP", "1"))
LEV = int(os.environ.get("LEV", "0"))
solver = os.environ.get("SOLVER", "nhflow")
n, nz, dx = 40, 6, 0.05
L = n*dx
h = n//2
I, J = np.meshgrid(np.arange(n+1), np.arange(n+1), indexing="ij")
x, y = I*dx, J*dx
X = x + 0.15*y + 0.04*np.sin(2*np.pi*x/L)*np.sin(np.pi*y/L)
Y = y + 0.05*x + 0.03*np.sin(np.pi*x/L)*np.sin(2*np.pi*y/L)
Z = float(os.environ.get("BUMP", "0.06"))*np.exp(-((X-1.3)**2+(Y-0.7)**2)/0.04)

def quad(i0, j0, k):
    return tuple(np.rot90(A[i0:i0+h+1, j0:j0+h+1], k=-k).copy() for A in (X, Y, Z))

def rot(m, k, a, b):
    """cell (a,b) of an m x m block turned by k quarter turns: its cell in the unturned block"""
    q = np.rot90(np.arange(m*m).reshape(m, m), k=-k)[a, b]
    return q//m, q % m

# cell areas of the global grid (shoelace)
def area(X, Y):
    x1, y1, x2, y2 = X[:-1, :-1], Y[:-1, :-1], X[1:, :-1], Y[1:, :-1]
    x3, y3, x4, y4 = X[1:, 1:], Y[1:, 1:], X[:-1, 1:], Y[:-1, 1:]
    return 0.5*np.abs((x3-x1)*(y4-y2) - (x4-x2)*(y3-y1))
AR = area(X, Y)

offs = [(0, 0, 0), (h, 0, 1), (0, h, 2), (h, h, 3)] if not os.environ.get("NOTURN") else [(0, 0, 0), (h, 0, 0), (0, h, 0), (h, h, 0)]
blocks = [quad(*o) for o in offs]
control = "B 1 %g\nB 10 0.0 %g 0.0 %g 0.0 %g\n\nM 10 NP\n\nR 20 1\n" % (dx, L, L, nz*dx)
# SFLOW: A 209 0 (the bed smoothing of SFLOW acts on the native grid only, not on the blocks)
flow = ("A 10 5\nA 520 %s\n" % os.environ.get("A520", "0")) if solver == "nhflow" else "A 10 2\nA 209 0\nA 210 3\nA 211 4\nA 220 0\nA 240 1\n"
sed = "S 10 1\nS 11 2\nS 20 0.0005\nS 15 3\nS 17 1.0\nS 35 20.0\nS 41 1\nS 43 1\nS 44 1\nS 31 1\n"
amr = "G 1 %d\nG 2 0\nG 3 2\nG 4 8\nG 10 0.7 1.4 0.7 1.4\n" % LEV if LEV > 0 else ""
ctrl = (flow + "F 60 0.2\nF 72 0.65 1.45 0.65 1.45 0.38\nG 50 1\nG 52 1\nA 509 0\nN 47 0.3\nN 48 0\nN 49 0.005\n"
        "M 10 NP\nP 10 0\nW 22 -9.81\n" + sed + amr + os.environ.get("XCTRL", "").replace(";", "\n") + "\n")
ROTONE = bool(os.environ.get("ROTONE"))

def grids(d, npr):
    """the AMR dump of the run; a run without mesh refinement (SFLOW, one block, no box): the native grid
    from the regression dump (eta, WL, bed; i outer)"""
    if os.path.exists(os.path.join(d, "amr", "amr_r0.bin")):
        return amr_grids2(d, npr)
    out = []
    for r in range(npr):
        st = read_state(sorted(glob.glob(os.path.join(d, "dump", "state_*_r%d.bin" % r)))[-1])
        G = {"rank": r, "id": 0, "lev": 0, "blk": 0, "I0": 0, "J0": 0, "nx": n, "ny": n, "count": st["count"], "simtime": st["simtime"]}
        for k in ("eta", "WL", "bed"):
            G[k] = np.asarray(st[k]).reshape(n, n)
        G["U"] = G["V"] = np.full((n, n), np.nan)
        out.append(G)
    return out

def cell_key(G4, g, a, b, kind):
    """the global cell (level, I, J) of cell (a,b) of grid g, None outside; kind one|four|rotone. The blocks 1..
    and the patches of the four-block run are in their block frame (offset: the level-0 box of the block)"""
    lv = g["lev"]
    m = h << lv
    if kind == "one":
        I0, J0 = g["I0"]+a, g["J0"]+b
        return (lv, I0, J0) if (I0 < n << lv and J0 < n << lv) else None
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

def table(G4, kind):
    t = {}
    for g in G4:
        for a in range(g["nx"]):
            for b in range(g["ny"]):
                key = cell_key(G4, g, a, b, kind)
                if key is not None:
                    t[key] = {f: g[f][a, b] for f in ("bed", "eta", "U", "V")}
    return t

RUNS = [("one0", [(X, Y, Z)], 1, 1, "one"), ("one", [(X, Y, Z)], 1, steps, "one"), ("four", blocks, np_, steps, "four")]
if ROTONE:
    RUNS[2] = ("four", [tuple(np.rot90(A, k=-1).copy() for A in (X, Y, Z))], np_, steps, "rotone")
T = {}
for name, bl, npr, st, kind in RUNS:
    d = os.path.join(S5, "runs/t5_%s_%s_L%d_np%d" % (solver, name, LEV, npr))
    c = ctrl + ("G 1 1\n" if len(bl) == 1 and LEV == 0 else "") + ("S 35 0.0\n" if name == "one0" else "")
    rc = 0
    if not os.environ.get("NORUN"):
        rc = run4(d, control.replace("NP", str(npr)), c.replace("NP", str(npr)), bl, np_=npr, steps=st)
    if rc != 0:
        print(name, "rc", rc, log_lines(d, "!!!")[:3])
    G = grids(d, npr)
    T[name] = table(G, kind)
    if name != "one0":
        print("%-4s rc=%d t=%.3f grids %d (levels %s)" % (name, rc, G[0]["simtime"], len(G), sorted({g["lev"] for g in G})))

def level0(t, f):
    a = np.full((n, n), np.nan)
    for (lv, I0, J0), v in t.items():
        if lv == 0:
            a[I0, J0] = v[f]
    return a

b0, b1, b4 = level0(T["one0"], "bed"), level0(T["one"], "bed"), level0(T["four"], "bed")
dv1, dv4 = (AR*(b1-b0)).sum(), (AR*(b4-b0)).sum()
ref = max((AR*np.abs(b1-b0)).sum(), 1e-300)
print("bed change: max %.3e; volume change one %.2e four %.2e (rel %.1e %.1e)" % (np.nanmax(np.abs(b1-b0)), dv1, dv4, abs(dv1)/ref, abs(dv4)/ref))
for lv in sorted({k[0] for k in T["one"]}):
    common = [k for k in T["one"] if k[0] == lv and k in T["four"]]
    nmiss = sum(1 for k in T["one"] if k[0] == lv) - len(common)
    err = {f: max([abs(T["four"][k][f]-T["one"][k][f]) for k in common if T["one"][k][f] == T["one"][k][f]] or [float("nan")]) for f in ("bed", "eta", "U", "V")}
    print("four vs one (level %d, %d cells%s): " % (lv, len(common), ", %d missing" % nmiss if nmiss else "")
          + "  ".join("%s %.1e" % (f, e) for f, e in err.items()))
if os.environ.get("TOP"):
    e = sorted(((abs(T["four"][k]["bed"]-T["one"][k]["bed"]), k) for k in T["one"] if k in T["four"]), reverse=True)[:12]
    print("largest |dbed| at", ["%s %.1e" % (k, v) for v, k in e])
