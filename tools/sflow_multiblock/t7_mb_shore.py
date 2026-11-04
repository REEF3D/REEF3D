# P2b T7: wetting and drying with refinement: dam break running up a slope across a rotated block
# interface, refinement of the shoreline (G 22) and of surface jumps (G 20), regrid every 4 steps:
# one block vs. two blocks, level 0 and the fine cells compared at the same cell centres
import sys, os, glob, numpy as np
from p2lib import *
from vtsfull import readvts_full
steps = int(sys.argv[1]) if len(sys.argv) > 1 else 60
np_ = int(sys.argv[2]) if len(sys.argv) > 2 else 1
lev = int(os.environ.get("LEV", "2"))
n, dx = 64, 0.5
I, J = np.meshgrid(np.arange(n+1), np.arange(n+1), indexing="ij")
X = I*dx + 0.1*J*dx; Y = J*dx + 0.3*np.sin(np.pi*I/n)*np.sin(2*np.pi*J/n)
Z = 0.08*Y + 0.3*np.exp(-((X-20.0)**2+(Y-10.0)**2)/8.0)
h = n//2
def rot_right(A):
    return np.array([[A[n-b, a] for b in range(h+1)] for a in range(n+1)])
two = [(X[:h+1], Y[:h+1], Z[:h+1]), (rot_right(X), rot_right(Y), rot_right(Z))]
ctrl = ("A 10 2\nA 210 3\nA 211 4\nA 220 0\nA 240 1\nA 209 0\n\nF 60 1.0\nF 72 10.0 22.0 2.0 8.0 2.0\n\n"
        "G 50 1\nG 1 %d\nG 2 4\nG 3 2\nG 4 8\nG 20 0.02\nG 22 1\n\nN 47 0.3\n\nM 10 %%d\n\nP 10 1\nP 181 %d\n\nW 22 -9.81\n" % (lev, steps))
control = "B 1 0.5\nB 10 0.0 32.0 0.0 32.0 0.0 0.5\n\nM 10 %d\nM 20 2\n\nR 20 1\n"
def cells(d):
    out = {}
    for f in sorted(glob.glob(os.path.join(d, "REEF3D_SFLOW_AMR", "REEF3D-SFLOW-AMR-00000001-*-*.vts"))):
        dd, xc, yc, lv = readvts_full(f)
        for a, b, e in zip(xc.ravel(), yc.ravel(), dd["eta"].ravel()):
            out[(lv, round(a, 7), round(b, 7))] = e
    return out
d1 = os.path.join(S2, "runs/t7_one"); d2 = os.path.join(S2, "runs/t7_two_np%d" % np_)
print("one rc", run_mb(d1, control % 1, ctrl % 1, [(X, Y, Z)], steps=steps))
print("two rc", run_mb(d2, control % np_, ctrl % np_, two, np_=np_, steps=steps), log_lines(d2, "!!!|nesting|source")[:3])
A = cells(d1); B = cells(d2)
fineA = {k[1:]: v for k, v in A.items() if k[0] >= 1}; fineB = {k[1:]: v for k, v in B.items() if k[0] >= 1}
c = set(fineA) & set(fineB)
print("fine cells %d / %d, common %d, max diff %.2e" % (len(fineA), len(fineB), len(c), max([abs(fineA[k]-fineB[k]) for k in c] or [0])))
ref = final_state(d1)["eta"].reshape(n, n)
if np_ == 1:
    s0 = final_state(d2)["eta"].reshape(h+8, n)
    print("block 0 level 0: max diff %.2e" % np.abs(s0[:h]-ref[:h]).max())
for d in (d1, d2):
    print(os.path.basename(d), open(os.path.join(d, "REEF3D_SFLOW_AMR", "REEF3D_SFLOW_AMR_log.dat")).read().splitlines()[-1])
s = final_state(d1); print("dry cells (one block):", int((s["WL"] < 1e-3).sum()))
