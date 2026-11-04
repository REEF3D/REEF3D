# P5a T4: sediment with mesh refinement (NHFLOW G 1): the bed-load dam break of T1 with a static refinement box
# over the column (LEV levels) or patches following the surface (ADAPT=1: G 20, regrid every 4 steps):
# Cartesian, curvilinear on the same nodes, turned by 90 deg; level 0 and every fine cell per cell index, the
# bed volume over the leaf cells (it starts flat: the sum of area times bed stays 0).
# environment: SOLVER nhflow|sflow, A520, LEV, ADAPT, XCTRL (;), NP
import sys, numpy as np
from p5lib import *
sys.path.insert(0, os.path.join(S5, "../nhflow_multiblock"))
from p4lib import amr_grids as _ag
steps = int(sys.argv[1]) if len(sys.argv) > 1 else 40
nx, ny, nz, dx = 32, 24, 6, 0.125
np_ = int(os.environ.get("NP", "1"))
LEV = int(os.environ.get("LEV", "1"))
ADAPT = int(os.environ.get("ADAPT", "0"))
ctl = control(nx, ny, nz, nx*dx, ny*dx, 0.5, np_=np_)
sed = "S 10 1\nS 11 2\nS 20 0.0005\nS 15 3\nS 17 1.0\nS 35 20.0\nS 41 1\nS 43 1\nS 44 1\nS 31 1\n"
amr = "G 1 %d\nG 3 2\nG 4 8\n" % LEV
box = (1.0, 1.75, 1.0, 2.0)
rbox = (0.5, 2.5, 0.75, 2.25)
X0, Y0 = 1000.0, 2000.0
solver = os.environ.get("SOLVER", "nhflow")
flow = ("A 10 5\nA 520 %s\n" % os.environ.get("A520", "0")) if solver == "nhflow" else "A 10 2\nA 210 3\nA 211 4\nA 220 0\nA 240 1\n"
base = (flow + "N 47 0.3\nN 48 0\nN 49 0.01\nF 60 0.2\nM 10 %d\nP 10 0\nW 22 -9.81\n" % (np_)
        + sed + amr + os.environ.get("XCTRL", "").replace(";", "\n"))
res = {}
for name, rot in [("cart", None), ("curv_id", "id"), ("rot90", "rot90")]:
    if rot is None:
        grid, b, rb = None, box, rbox
    else:
        X, Y = nodes(nx, ny, dx, rot, X0, Y0)
        grid = (X, Y, np.zeros_like(X))
        b, rb = box_phys(box, rot, X0, Y0), box_phys(rbox, rot, X0, Y0)
    ctrl = base + "F 72 %g %g %g %g 0.35\n" % b + ("G 2 4\nG 20 0.01\n" if ADAPT else "G 2 0\nG 10 %g %g %g %g\n" % rb) + ("G 50 2\n" if grid is not None else "")
    d = os.path.join(S5, "runs/t4_%s_%s_L%d_A%d" % (solver, name, LEV, ADAPT))
    rc = run5(d, ctl, ctrl, grid=grid, np_=np_, steps=steps)
    G = amr_grids2(d, np_)
    res[name] = G
    # leaf volume: cells of level l not covered by level l+1
    cov = set()
    for g in G:
        if g["lev"] >= 1:
            for a in range(g["nx"]):
                for bb in range(g["ny"]):
                    cov.add((g["lev"]-1, (g["I0"]+a)//2, (g["J0"]+bb)//2))
    vol, n = 0.0, {}
    bmax = 0.0
    for g in G:
        l = g["lev"]
        A = (dx/2**l)**2
        n[l] = n.get(l, 0) + g["nx"]*g["ny"]
        for a in range(g["nx"]):
            for bb in range(g["ny"]):
                if (l, g["I0"]+a, g["J0"]+bb) not in cov:
                    vol += A*g["bed"][a, bb]
        bmax = max(bmax, np.abs(g["bed"]).max())
    print("%-8s rc=%d t=%.3f  cells per level %s  max|bed| %.3e  leaf volume %.2e (rel %.1e)" % (name, rc, G[0]["simtime"], dict(sorted(n.items())), bmax, vol,
          abs(vol)/max(bmax*nx*ny*dx*dx, 1e-30)), log_lines(d, "!!!")[:2])
def table(G):
    t = {}
    for g in G:
        for a in range(g["nx"]):
            for bb in range(g["ny"]):
                t[(g["lev"], g["I0"]+a, g["J0"]+bb)] = (g["bed"][a, bb], g["eta"][a, bb])
    return t
T = {k: table(v) for k, v in res.items()}
for name, ref in [("curv_id", "cart"), ("rot90", "curv_id")]:
    a, c = T[name], T[ref]
    if set(a) != set(c):
        print("%-8s vs %-7s: different cells" % (name, ref))
        continue
    eb = max(abs(a[k][0]-c[k][0]) for k in a)
    ee = max(abs(a[k][1]-c[k][1]) for k in a)
    print("%-8s vs %-7s: bed %.1e  eta %.1e  (all levels)" % (name, ref, eb, ee))
