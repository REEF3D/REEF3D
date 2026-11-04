# P2 T1: dam break in a closed basin across a block interface: one block vs. two blocks, the second
# with rotated indices (its i along +y, j along -x): the same eta per cell
import sys, numpy as np
from p2lib import *
from vtsfull import readvts
import glob

steps = int(sys.argv[1]) if len(sys.argv) > 1 else 40
np_ = int(sys.argv[2]) if len(sys.argv) > 2 else 1
n, dx = 60, 0.5
I, J = np.meshgrid(np.arange(n+1), np.arange(n+1), indexing="ij")
X, Y = I*dx, J*dx
Z = 0.3*np.exp(-((X-20.0)**2+(Y-10.0)**2)/8.0)
ctrl = "A 10 2\nA 210 3\nA 211 4\nA 220 0\nA 240 1\nA 209 0\n\nF 60 1.0\nF 72 12.0 18.0 12.0 18.0 2.0\n\nG 50 1\n\nN 47 0.3\n\nM 10 %d\n\nP 10 1\nP 181 %d\n\nW 22 -9.81\n" % (np_, steps)
control = "B 1 0.5\nB 10 0.0 30.0 0.0 30.0 0.0 0.5\n\nM 10 %d\nM 20 2\n\nR 20 1\n" % np_

# block 1: index (a,b) -> node (60-b, a) of the right half
def rot_right(A):
    return np.array([[A[60-b, a] for b in range(31)] for a in range(61)])
one = [(X, Y, Z)]
two = [(X[:31], Y[:31], Z[:31]), (rot_right(X), rot_right(Y), rot_right(Z))]

R = {}
for name, blocks in [("one", one), ("two", two)]:
    npr = 1 if name == "one" else np_
    d = os.path.join(S2, "runs/t1_%s_np%d" % (name, npr))
    rc = run_mb(d, control.replace("M 10 %d" % np_, "M 10 %d" % npr), ctrl.replace("M 10 %d" % np_, "M 10 %d" % npr), blocks, np_=npr, steps=steps)
    print(name, "rc", rc, [l for l in log_lines(d, "!!!|multi-block")][:4] if rc >= -1 else open(os.path.join(d,"dive.log")).read()[-500:])
    R[name] = d

ref = final_state(R["one"])
eref = ref["eta"].reshape(n, n)
st = [final_state(R["two"], r) for r in range(np_)]
# native (block 0 with its ring, 38 x 60) from the regression dump, block 1 from the AMR VTS output
e0 = np.concatenate([s["eta"] for s in st]) if np_ == 1 else None
vts = sorted(glob.glob(os.path.join(R["two"], "REEF3D_SFLOW_AMR", "REEF3D-SFLOW-AMR-[0-9]*-0001-0001.vts")))
d1, nx1, ny1 = readvts(vts[-1])
e1 = d1["eta"]
diff1 = max(abs(e1[a, b] - eref[59-b, a]) for a in range(nx1) for b in range(ny1))
print("block 1 (from VTS, %s): max |eta - one block| = %.2e" % (os.path.basename(vts[-1]), diff1))
if e0 is not None:
    e0 = e0.reshape(38, n)
    print("block 0: max |eta - one block| = %.2e, ring %.2e" % (np.abs(e0[:30]-eref[:30]).max(), np.abs(e0[30:]-eref[30:38]).max()))
print("eta range one block %.3f .. %.3f" % (eref.min(), eref.max()))
