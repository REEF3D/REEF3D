# P4 T1: NHFLOW dam break in a closed basin across a block interface: one block vs. two blocks, the
# second with its indices turned (its i along +y, j along -x): per cell the same state (velocities are
# Cartesian components); A520="0 2" the pressure schemes, NP ranks for the two-block run
import sys, numpy as np
from p4lib import *

steps = int(sys.argv[1]) if len(sys.argv) > 1 else 30
np_ = int(sys.argv[2]) if len(sys.argv) > 2 else 1
schemes = os.environ.get("A520", "0 2").split()
n, nz, dx = 40, 8, 0.05
h = n//2
I, J = np.meshgrid(np.arange(n+1), np.arange(n+1), indexing="ij")
X, Y = I*dx, J*dx
Z = float(os.environ.get("BUMP","-0.1"))*np.exp(-((X-1.3)**2+(Y-0.6)**2)/0.05)
XCTRL = os.environ.get("XCTRL", "").replace(";", "\n")
control = "B 1 %g\nB 10 0.0 %g 0.0 %g 0.0 %g\n\nM 10 NP\n\nR 20 1\n" % (dx, n*dx, n*dx, nz*dx)

def rot(A):     # block 1: index (a,b) -> node (n-b, a) of the right half
    return np.array([[A[n-b, a] for b in range(h+1)] for a in range(n+1)])
one = [(X, Y, Z)]
ROT = int(os.environ.get("ROT", "1"))
two = [(X[:h+1], Y[:h+1], Z[:h+1]), (rot(X), rot(Y), rot(Z)) if ROT else (X[h:], Y[h:], Z[h:])]

for a520 in schemes:
    ctrl = "A 10 5\nA 520 %s\nF 60 0.5\nF 72 0.91 1.11 0.71 0.91 0.56\nG 50 1\nG 52 1\nN 47 0.3\nM 10 NP\nP 10 0\nW 22 -9.81\n%s\n" % (a520, XCTRL)
    R = {}
    for name, blocks, npr in [("one", one, 1), ("two", two, np_)]:
        d = os.path.join(S4, "runs/t1_%s_A%s_np%d" % (name, a520, npr))
        # the one-block reference: the hierarchy without patches (G 1 1, nothing to refine) dumps its grid
        c = ctrl + ("G 1 1\n" if name == "one" else "")
        rc = run4(d, control.replace("NP", str(npr)), c.replace("NP", str(npr)), blocks, np_=npr, steps=steps)
        print("A520 %s %-3s rc %d" % (a520, name, rc), log_lines(d, "!!!|multi-block")[:3] if rc > -100 else open(os.path.join(d, "dive.log")).read()[-400:])
        R[name] = amr_grids(d, npr)
    ref = level0(R["one"], n, n)
    rg = R["two"]
    nat = level0(rg, h+8, n)
    blk1 = [G for G in rg if G["blk"] == 1][0]
    out = []
    for f in FIELDS:
        e0 = np.abs(nat[f][:h]-ref[f][:h]).max()
        er = np.abs(nat[f][h:h+8]-ref[f][h:h+8]).max()
        b = blk1[f]
        if ROT:
            e1 = max(np.abs(b[a, bb]-ref[f][n-1-bb, a]).max() for a in range(blk1["nx"]) for bb in range(blk1["ny"]))
        else:
            e1 = np.abs(b-ref[f][h:]).max()
        out.append("%s %.1e/%.1e/%.1e" % (f, e0, e1, er))
    print("A520 %s step %d t %.4f |eta| %.3f |U| %.3f |P| %.1f   block 0 / block 1 / ring: %s" % (a520, blk1["count"], blk1["simtime"], np.abs(ref["eta"]).max(), np.abs(ref["U"]).max(), np.abs(ref["P"]).max(), "  ".join(out)))
