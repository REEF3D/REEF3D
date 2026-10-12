# P2b T6: three blocks around a singular vertex (a hexagon split into three quadrilaterals), dam break
# over the vertex, refinement box over it: no refinement at the vertex, mass at round-off, the level-0
# surface close to the run without refinement, no spikes
import sys, os, glob, numpy as np
from p2lib import *
from vtsfull import readvts_full

steps = int(sys.argv[1]) if len(sys.argv) > 1 else 60
np_ = int(sys.argv[2]) if len(sys.argv) > 2 else 1
n = 32
R = 16.0
O = np.array([20.0, 20.0])
V = [O + R*np.array([np.cos(np.pi/3*k), np.sin(np.pi/3*k)]) for k in range(6)]
a = np.linspace(0, 1, n+1)
A, C = np.meshgrid(a, a, indexing="ij")
blocks = []
for b in range(3):
    P0, P1, P2, P3 = O, V[2*b], V[(2*b+1) % 6], V[(2*b+2) % 6]
    # node (i,j): i from O towards V[2b], j from O towards V[2b+2] (right-handed)
    X = (1-A)*(1-C)*P0[0] + A*(1-C)*P1[0] + (1-A)*C*P3[0] + A*C*P2[0]
    Y = (1-A)*(1-C)*P0[1] + A*(1-C)*P1[1] + (1-A)*C*P3[1] + A*C*P2[1]
    Z = 0.2*np.exp(-((X-26.0)**2+(Y-22.0)**2)/10.0)
    blocks.append((X, Y, Z))
ctrl = ("A 10 2\nA 210 3\nA 211 4\nA 220 0\nA 240 1\nA 209 0\n\nF 60 1.0\nF 72 17.0 23.0 17.0 23.0 1.6\n\n"
        "G 50 1\nG 1 %d\nG 2 4\nG 3 2\nG 4 8\nG 10 10.0 30.0 10.0 30.0\nG 20 0.02\n\nN 47 0.3\n\nM 10 %%d\n\nP 10 1\nP 181 %d\n\nW 22 -9.81\n")
control = "B 1 0.5\nB 10 0.0 16.0 0.0 16.0 0.0 0.5\n\nM 10 %d\nM 20 2\n\nR 20 1\n"
res = {}
for lev in (0, 2):
    d = os.path.join(S2, "runs/t6_L%d_np%d" % (lev, np_))
    rc = run_mb(d, control % np_, ctrl % (lev, steps) % np_, blocks, np_=np_, steps=steps)
    print("levels", lev, "rc", rc, log_lines(d, "!!!|SFLOW AMR:|singular|without a source|nesting")[:5])
    lg = open(os.path.join(d, "REEF3D_SFLOW_AMR", "REEF3D_SFLOW_AMR_log.dat")).read().splitlines()
    print("   last log:", lg[-1])
    cells0, fine = {}, []
    for f in sorted(glob.glob(os.path.join(d, "REEF3D_SFLOW_AMR", "REEF3D-SFLOW-AMR-00000001-*-*.vts"))):
        dd, xc, yc, lv = readvts_full(f)
        if lv == 0:
            for x_, y_, e in zip(xc.ravel(), yc.ravel(), dd["eta"].ravel()):
                cells0[(round(x_, 6), round(y_, 6))] = e
        else:
            r = np.sqrt((xc-O[0])**2 + (yc-O[1])**2)
            fine.append((lv, r.min(), np.nanmax(np.abs(dd["eta"])), np.isnan(dd["eta"]).any()))
    if np_ == 1:
        st = final_state(d)
        res[lev] = (cells0, st["eta"])
    if fine:
        print("   patches: %d, closest fine cell to the vertex %.2f m (level-0 spacing ~%.2f m), max |eta| %.3f, NaN %s" %
              (len(fine), min(f[1] for f in fine), R/n, max(f[2] for f in fine), any(f[3] for f in fine)))
if np_ == 1:
    c0, c2 = res[0][0], res[2][0]
    k = set(c0) & set(c2)
    print("blocks 1,2 level 0: max |eta(G 1 2) - eta(no refinement)| = %.2e over %d cells" % (max(abs(c0[q]-c2[q]) for q in k), len(k)))
    print("block 0: max |diff| = %.2e" % np.abs(res[0][1]-res[2][1]).max())
