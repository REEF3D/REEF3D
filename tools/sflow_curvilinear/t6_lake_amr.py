# T6: lake at rest on the wavy, skewed grid with patches: static box patches (G 10) and patches that
# are rebuilt every 4 steps (G 2 4, flagged by a box and by the shoreline G 22), hydrostatic,
# non-hydrostatic and Boussinesq (one level: dx/h >= 1/4 on the skewed grid), with and without
# subcycling; velocities must stay at round-off,
# the water volume constant (area-weighted restriction and prolongation)
import sys, numpy as np
from p1lib import *

steps = int(sys.argv[1]) if len(sys.argv) > 1 else 80
n, dx = 60, 0.5
L = n*dx
control = "B 1 %g\nB 10 0.0 %g 0.0 %g 0.0 0.5\n\nM 10 1\nM 20 2\n" % (dx, L, L)
I, J = np.meshgrid(np.arange(n+1), np.arange(n+1), indexing="ij")
x, y = I*dx, J*dx
X = x + 0.25*y + 0.8*np.sin(2*np.pi*x/L)*np.sin(np.pi*y/L)
Y = y + 0.1*x + 0.6*np.sin(np.pi*x/L)*np.sin(4*np.pi*y/L)
rng = np.random.default_rng(1)
Zsub = 0.6*np.exp(-((X-17.0)**2+(Y-17.5)**2)/12.0) + 0.05*rng.standard_normal(X.shape)

base = "A 10 2\nA 210 3\nA 211 4\nA 240 1\nF 60 1.0\n\nG 50 2\n\nN 41 100.0\nN 47 0.3\n\nM 10 1\n\nW 22 -9.81\n\nG 3 2\nG 4 8\n%s%s"
amr = [("static", "G 2 0\nG 10 6.0 20.0 8.0 22.0\n"), ("regrid", "G 2 4\nG 10 6.0 20.0 8.0 22.0\nG 22 1\n")]
eqs = [("hyd", "A 220 0\nG 1 2\n"), ("nh", "A 220 2\nG 1 2\n"), ("bous", "A 220 4\nG 1 1\n"), ("nh_sub", "A 220 2\nG 1 2\nG 7 1\n"), ("bous_sub", "A 220 4\nG 1 1\nG 7 1\n")]

sel = sys.argv[2:]
for an, ak in amr:
    for en, ek in eqs:
        if sel and an+"_"+en not in sel:
            continue
        dd = os.path.join(S, "runs/t6_%s_%s" % (an, en))
        rc = run(dd, control, base % (ek, ak), grid=(X, Y, Zsub), steps=steps)
        if rc != 0:
            print(an, en, "FAILED", log_lines(dd, "!!!")[:3]); continue
        st = final_state(dd)
        rows = [l.split() for l in open(os.path.join(dd, "REEF3D_SFLOW_AMR", "REEF3D_SFLOW_AMR_log.dat")) if l.strip() and l.split()[0].isdigit()]
        print("%-7s %-7s patches %s cells %s | max |P| %.2e |Q| %.2e |eta| %.2e | mass rel. change %s"
              % (an, en, rows[-1][3], rows[-1][4], np.abs(st["P"]).max(), np.abs(st["Q"]).max(), np.abs(st["eta"]).max(), rows[-1][6]))
