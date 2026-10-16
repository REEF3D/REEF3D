# P3 T8: mesh refinement on a wavy, skewed curvilinear grid (non-affine cells): (a) a lake at rest over a
# bed bump with a static patch across the bump: all velocities and the pressure at round-off on both
# levels; (b) a dam break with patches following the surface (regrid every step, G 7 0 / 1): the water
# volume of the composite grid to round-off (conservative restriction, refluxing per face length)
import sys, glob, numpy as np
from p3lib import *
import p1lib
from readvts import readvts
steps = int(sys.argv[1]) if len(sys.argv) > 1 else 40
n, nz, dx = 48, 6, 0.5
L = n*dx
control = "C 11 21\nC 12 21\nC 13 21\nC 14 21\nC 15 21\nC 16 3\n\nB 2 %d %d %d\nB 10 0.0 %g 0.0 %g 0.0 1.0\n\nM 10 1\n" % (n, n, nz, L, L)
I, J = np.meshgrid(np.arange(n+1), np.arange(n+1), indexing="ij")
x, y = I*dx, J*dx
if os.environ.get("FLAT","0")=="1":
    X, Y = x.astype(float), y.astype(float)
else:
  X = x + 0.25*y + 0.8*np.sin(2*np.pi*x/L)*np.sin(np.pi*y/L)
  Y = y + 0.1*x + 0.6*np.sin(np.pi*x/L)*np.sin(4*np.pi*y/L)
def run(name, ctrl, Z):
    d = os.path.join(S3, "runs/t8_"+name)
    rc = p1lib.run(d, control, ctrl, grid=(X, Y, Z), steps=steps)
    st = final_state(d)
    log = [l.split() for l in open(os.path.join(d, "REEF3D_NHFLOW_AMR", "REEF3D_NHFLOW_AMR_log.dat")) if l.strip() and l.split()[0].isdigit()]
    last = sorted(glob.glob(os.path.join(d, "REEF3D_NHFLOW_AMR", "REEF3D-NHFLOW-AMR-*.vtm")))[-1]
    stamp = os.path.basename(last).split("-")[-1].split(".")[0]
    pf = sorted(f for f in glob.glob(os.path.join(d, "REEF3D_NHFLOW_AMR", "REEF3D-NHFLOW-AMR-%s-*" % stamp)) if not f.endswith(".vtm") and "-L0-" not in f)
    pats = [readvts(f)[0] for f in pf]
    msg = [l for l in p1lib.log_lines(d, "!!!")]
    return rc, st, log, pats, msg
# (a) lake at rest, bed bump under the patch, static box over the middle (physical coordinates)
Zb = 0.6*np.exp(-((X-14.0)**2+(Y-13.0)**2)/10.0)
for A in [0, 2]:
    c = "A 10 5\nA 520 %d\nF 60 1.0\nN 47 0.3\nM 10 1\nP 10 0\nP 30 0.5\nW 22 -9.81\nG 50 2\nG 1 1\nG 10 9.0 19.0 8.0 18.0\nG 4 4\n" % A
    rc, st, log, pats, msg = run("lake_A%d" % A, c, Zb)
    pu = max([np.abs(p["u_top"]).max() for p in pats] + [np.abs(p["v_top"]).max() for p in pats] + [0.0])
    pe = max([np.abs(p["eta"]).max() for p in pats] + [0.0])
    print("lake A 520 %d rc=%d patches %s: level 0 max |U| %.1e |V| %.1e |W| %.1e |P| %.1e |eta| %.1e; patches |u|,|v| %.1e |eta| %.1e" % (A, rc, log[-1][3], np.abs(st["U"]).max(), np.abs(st["V"]).max(), np.abs(st["W"]).max(), np.abs(st["P"]).max(), np.abs(st["eta"]).max(), pu, pe), msg)
# (b) dam break, adaptive, mass
Z0 = np.zeros_like(X)
for sub in [0, 1]:
    c = "A 10 5\nA 520 2\nF 60 1.0\nF 72 10.0 14.0 10.0 14.0 1.3\nN 47 0.3\nM 10 1\nP 10 0\nP 30 0.5\nW 22 -9.81\nG 50 2\nG 1 1\nG 20 0.005\nG 4 4\nG 2 1\nG 5 0\nG 7 %d\n" % sub
    rc, st, log, pats, msg = run("dam_sub%d" % sub, c, Z0)
    rel = [float(l[8]) for l in log]
    print("dam break G 7 %d rc=%d steps %s patches %s columns %s: max |volume change| %.1e" % (sub, rc, log[-1][0], log[-1][3], log[-1][4], np.abs(rel).max()), msg)
