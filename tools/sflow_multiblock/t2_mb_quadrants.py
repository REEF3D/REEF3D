# P2 T2: four blocks meeting at a regular vertex, each with another orientation of its indices
# (environment: SKEW 0 a Cartesian grid, RW the ring width R 21, XCTRL extra ctrl lines separated by ;)
# (0, 90, 180, 270 deg), dam break over the vertex, wavy skewed grid, bumpy bed: one block vs. four
import sys, glob, numpy as np
from p2lib import *
from vtsfull import readvts

steps = int(sys.argv[1]) if len(sys.argv) > 1 else 40
np_ = int(sys.argv[2]) if len(sys.argv) > 2 else 1
n, dx = 60, 0.5
L = n*dx
I, J = np.meshgrid(np.arange(n+1), np.arange(n+1), indexing="ij")
x, y = I*dx, J*dx
SK = float(os.environ.get("SKEW","1"))
X = x + SK*0.15*y + SK*0.5*np.sin(2*np.pi*x/L)*np.sin(np.pi*y/L)
Y = y + SK*0.05*x + SK*0.4*np.sin(np.pi*x/L)*np.sin(2*np.pi*y/L)
Z = 0.3*np.exp(-((X-20.0)**2+(Y-12.0)**2)/8.0)
h = n//2

# a block of the quadrant (i0..i0+h, j0..j0+h) with its indices turned by k*90 deg:
# node (a,b) of the block -> global node
def quad(i0, j0, k):
    out = []
    for A in (X, Y, Z):
        Q = A[i0:i0+h+1, j0:j0+h+1]          # (h+1, h+1) index [di, dj]
        out.append(np.rot90(Q, k=-k).copy())  # k quarter turns of the index frame
    return tuple(out)

def gidx(i0, j0, k, a, b):
    """global cell of cell (a,b) of the block quad(i0,j0,k)"""
    M = np.arange(h*h).reshape(h, h)
    R = np.rot90(M, k=-k)
    q = R[a, b]
    return i0 + q//h, j0 + q % h

blocks = [quad(0, 0, 0), quad(h, 0, 1), quad(0, h, 2), quad(h, h, 3)]
offs = [(0, 0, 0), (h, 0, 1), (0, h, 2), (h, h, 3)]
ctrl = "A 10 2\nA 210 3\nA 211 4\nA 220 0\nA 240 1\nA 209 0\n\nF 60 1.0\nF 72 12.0 18.0 12.0 18.0 2.0\n\nG 50 1\n\nN 47 0.3\n\nM 10 %d\n\nP 10 1\nP 181 %d\n\nW 22 -9.81\n" + os.environ.get("XCTRL","").replace(";","\n") + "\n"
RW = int(os.environ.get("RW","8"))
control = "B 1 0.5\nB 10 0.0 30.0 0.0 30.0 0.0 0.5\n\nM 10 %d\nM 20 2\n\nR 20 1\nR 21 " + str(RW) + "\n"

d1 = os.path.join(S2, "runs/t2_one")
print("one rc", run_mb(d1, control % 1, ctrl % (1, steps), [(X, Y, Z)], steps=steps))
eref = final_state(d1)["eta"].reshape(n, n)
d4 = os.path.join(S2, "runs/t2_four_np%d" % np_)
rc = run_mb(d4, control % np_, ctrl % (np_, steps), blocks, np_=np_, steps=steps)
print("four rc", rc, log_lines(d4, "!!!|multi-block|without a source"))
print(open(os.path.join(d4, "DIVEMesh_Log/multiblock.log")).read())
vts = sorted(glob.glob(os.path.join(d4, "REEF3D_SFLOW_AMR", "REEF3D-SFLOW-AMR-00000001-*-*.vts")))
for f in vts:
    dd, nx, ny = readvts(f)
    # the block of the file: the orientation with the smallest difference
    best = None
    for bq in (1, 2, 3):
        i0, j0, k = offs[bq]
        err = max(abs(dd["eta"][a, b] - eref[gidx(i0, j0, k, a, b)]) for a in range(nx) for b in range(ny))
        best = err if best is None else min(best, err)
    print(os.path.basename(f), "max |eta - one block| = %.2e" % best)
if np_ == 1:
    e0 = final_state(d4)["eta"].reshape(h+RW, h+RW)
    print("block 0: %.2e" % np.abs(e0[:h, :h] - eref[:h, :h]).max())
print("mass", open(os.path.join(d4, "REEF3D_SFLOW_AMR", "REEF3D_SFLOW_AMR_log.dat")).read().splitlines()[-1].split()[-1])
