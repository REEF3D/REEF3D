# P2b T5: dynamic refinement (G 20 surface jumps, regrid every 4 steps) moving across the interfaces: one block vs. four
# (orientations 0, 90, 180, 270 deg), wavy skewed grid, bumpy bed, dam break over the vertex, a static
# refinement box over the vertex (all four interfaces): one block with the patches vs. four blocks
# with the patches; level 0 and every fine cell compared at the same cell centres
import sys, os, glob, numpy as np
from p2lib import *
from vtsfull import readvts, readvts_full

steps = int(sys.argv[1]) if len(sys.argv) > 1 else 20
np_ = int(sys.argv[2]) if len(sys.argv) > 2 else 1
lev = int(os.environ.get("LEV", "1"))
place = int(os.environ.get("PLACE", "0"))
n, dx = 64, 0.5
L = n*dx
SK = float(os.environ.get("SKEW", "1"))
I, J = np.meshgrid(np.arange(n+1), np.arange(n+1), indexing="ij")
x, y = I*dx, J*dx
X = x + SK*0.15*y + SK*0.5*np.sin(2*np.pi*x/L)*np.sin(np.pi*y/L)
Y = y + SK*0.05*x + SK*0.4*np.sin(np.pi*x/L)*np.sin(2*np.pi*y/L)
Z = 0.3*np.exp(-((X-20.0)**2+(Y-12.0)**2)/8.0)
h = n//2

def quad(i0, j0, k):
    return tuple(np.rot90(A[i0:i0+h+1, j0:j0+h+1], k=-k).copy() for A in (X, Y, Z))

blocks = [quad(0, 0, 0), quad(h, 0, 1), quad(0, h, 2), quad(h, h, 3)]
xc0 = 0.25*(X[:-1,:-1]+X[1:,:-1]+X[:-1,1:]+X[1:,1:]); yc0 = 0.25*(Y[:-1,:-1]+Y[1:,:-1]+Y[:-1,1:]+Y[1:,1:])
box = os.environ.get("TOL", "0.02")
ctrl = ("A 10 2\nA 210 3\nA 211 4\nA 220 0\nA 240 1\nA 209 0\n\nF 60 1.0\nF 72 14.0 20.0 13.0 19.0 2.0\n\n"
        "G 50 1\nG 1 %d\nG 2 4\nG 3 2\nG 4 8\nG 20 %s\nG 40 %d\n\nN 47 0.3\n\nM 10 %%d\n\nP 10 1\nP 181 %d\n\nW 22 -9.81\n" % (lev, box, place, steps))
control = "B 1 0.5\nB 10 0.0 %g 0.0 %g 0.0 0.5\n\nM 10 %%d\nM 20 2\n\nR 20 1\n" % (L, L)

def cells(d, nproc):
    """all fine cells of the AMR output (level >= 1): {(x,y) rounded: eta}, and the level-0 blocks"""
    fine, blk = {}, []
    for f in sorted(glob.glob(os.path.join(d, "REEF3D_SFLOW_AMR", "REEF3D-SFLOW-AMR-00000001-*-*.vts"))):
        dd, xc, yc, lv = readvts_full(f)
        if lv == 0:
            blk.append((dd["eta"], xc, yc))
            continue
        for a, b, e in zip(xc.ravel(), yc.ravel(), dd["eta"].ravel()):
            fine[(round(a, 7), round(b, 7))] = e
    return fine, blk

d1 = os.path.join(S2, "runs/t5_one")
print("one rc", run_mb(d1, control % 1, ctrl % 1, [(X, Y, Z)], steps=steps), log_lines(d1, "!!!|SFLOW AMR:")[:3])
d4 = os.path.join(S2, "runs/t5_four_np%d" % np_)
rc = run_mb(d4, control % np_, ctrl % np_, blocks, np_=np_, steps=steps)
print("four rc", rc, log_lines(d4, "!!!|SFLOW AMR:|multi-block|without a source|nesting")[:6])
ref = final_state(d1)["eta"].reshape(n, n)
fr, _ = cells(d1, 1)
f4, b4 = cells(d4, np_)
# level 0: blocks 1..3 by cell centre against the single-block level 0
key0 = {(round(a, 7), round(b, 7)): (i, j) for (i, j), a, b in zip(np.ndindex(n, n), xc0.ravel(), yc0.ravel())}
e0 = 0.0; miss = 0
for eta, xc, yc in b4:
    for a, b, e in zip(xc.ravel(), yc.ravel(), eta.ravel()):
        k = key0.get((round(a, 7), round(b, 7)))
        if k is None: miss += 1; continue
        e0 = max(e0, abs(e - ref[k]))
if np_ == 1:
    s0 = final_state(d4)["eta"].reshape(h+8, h+8)
    e0 = max(e0, np.abs(s0[:h, :h] - ref[:h, :h]).max())
print("level 0: max |eta - one block| = %.2e (cells not found %d)" % (e0, miss))
common = set(fr) & set(f4)
ef = max([abs(fr[k] - f4[k]) for k in common] or [0.0])
print("fine cells: one block %d, four blocks %d, common %d, max |eta diff| = %.2e" % (len(fr), len(f4), len(common), ef))
def lastmass(d):
    return open(os.path.join(d, "REEF3D_SFLOW_AMR", "REEF3D_SFLOW_AMR_log.dat")).read().splitlines()[-1].split()
print("one  log:", " ".join(lastmass(d1)))
print("four log:", " ".join(lastmass(d4)))
