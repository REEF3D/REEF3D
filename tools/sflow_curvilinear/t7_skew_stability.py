# T7: Boussinesq (A 220 4) on skewed curvilinear grids: a disturbance of 1e-6 in a basin at rest,
# flat bed, depth 1 m.  One level with the grid sheared up to 45 deg non-orthogonality, and with a
# refinement patch (finest dx = h/4) at 20 deg: sum(P^2) must stay bounded (it grew exponentially
# when the non-orthogonal terms of the own component were lagged in time)
import sys, numpy as np
import p1lib
from p1lib import S, run, log_lines

nsteps = int(sys.argv[1]) if len(sys.argv) > 1 else 60
base = "A 10 2\nA 210 3\nA 211 4\nA 220 4\nA 240 1\nF 60 1.0\nF 72 %s 1.000001\n\nG 50 2\n\nN 41 100.0\nN 47 0.3\n\nM 10 1\n\nW 22 -9.81\n%s"
amr = "G 1 1\nG 3 2\nG 4 8\nG 2 0\nG 10 %s\n"

def grid(n, dx, skew):
    I, J = np.meshgrid(np.arange(n+1), np.arange(n+1), indexing="ij")
    x, y = I*dx, J*dx
    return x + skew*0.25*y, y + skew*0.1*x

cases = [("one level, 20 deg", 1.0, False), ("one level, 32 deg", 2.0, False),
         ("one level, 45 deg", 2.5, False), ("patch dx=h/4, 20 deg", 1.0, True)]
n, dx = 60, 0.5
control = "B 1 %g\nB 10 0.0 %g 0.0 %g 0.0 0.5\n\nM 10 1\nM 20 2\n" % (dx, n*dx, n*dx)
for k, (name, skew, refine) in enumerate(cases):
    X, Y = grid(n, dx, skew)
    box = lambda a, b: "%g %g %g %g" % (X[a, a], X[b, b], Y[a, a], Y[b, b])
    c = base % (box(28, 32), (amr % box(20, 40)) if refine else "")
    dd = S + "/runs/t7_%d" % k
    rc = run(dd, control, c, grid=(X, Y, np.zeros_like(X)), steps=nsteps)
    if rc != 0:
        print(name, "FAILED", log_lines(dd, "!!!")[:3])
        continue
    a = p1lib.steps(dd)
    print("%-22s sum(P^2) step 10 %.2e, last %.2e, max %.2e" % (name, a[9, 3], a[-1, 3], a[:, 3].max()))
