# P3 T2: NHFLOW lake at rest with an island (wetting and drying) on a wavy, skewed curvilinear grid:
# velocities and non-hydrostatic pressure stay at round-off (A 520 0, 1, 2)
import sys, numpy as np
from p3lib import *
import p1lib
steps = int(sys.argv[1]) if len(sys.argv) > 1 else 100
n, nz, dx = 60, 6, 0.5
L = n*dx
control = "C 11 21\nC 12 21\nC 13 21\nC 14 21\nC 15 21\nC 16 3\n\nB 2 %d %d %d\nB 10 0.0 %g 0.0 %g 0.0 1.0\n\nM 10 1\n" % (n, n, nz, L, L)
I, J = np.meshgrid(np.arange(n+1), np.arange(n+1), indexing="ij")
x, y = I*dx, J*dx
X = x + 0.25*y + 0.8*np.sin(2*np.pi*x/L)*np.sin(np.pi*y/L)
Y = y + 0.1*x + 0.6*np.sin(np.pi*x/L)*np.sin(4*np.pi*y/L)
rng = np.random.default_rng(1)
Z = 1.3*np.exp(-((X-17.0)**2+(Y-17.5)**2)/12.0) + 0.05*rng.standard_normal(X.shape)
for A in [int(a) for a in os.environ.get("A520","0 2 1").split()]:
    ctrl = "A 10 5\nA 520 %d\nF 60 1.0\nN 47 0.3\nM 10 1\nP 10 0\nW 22 -9.81\nG 50 2\n" % A
    d = os.path.join(S3, "runs/t2_A%d" % A)
    rc = p1lib.run(d, control, ctrl, grid=(X, Y, Z), steps=steps)
    st = final_state(d)
    wl = st["WL"]; wet = wl > 0.002
    print("A 520 %d rc=%d count=%d t=%.3f wet columns %d of %d" % (A, rc, st["count"], st["simtime"], wet.sum(), wet.size))
    print("  max |U| %.2e |V| %.2e |W| %.2e |P| %.2e  max |eta| (wet) %.2e" % (np.abs(st["U"]).max(), np.abs(st["V"]).max(), np.abs(st["W"]).max(), np.abs(st["P"]).max(), np.abs(st["eta"][wet]).max()))
    for l in p1lib.log_lines(d, "!!!|nan|NaN"):
        print("  ", l)
