# P4 T3: the composite pressure over 1, 4, 9, 16 blocks (NB x NB, the indices of every block turned by
# another quarter turn): the same solution as one block, the pressure iterations against the block count
# (environment: NBS the block counts per direction, A520, G51 the passes of the block-local V-cycles)
import sys, numpy as np
from p4lib import *

steps = int(sys.argv[1]) if len(sys.argv) > 1 else 10
np_ = int(sys.argv[2]) if len(sys.argv) > 2 else 1
nbs = [int(v) for v in os.environ.get("NBS", "1 2 3 4").split()]
a520 = os.environ.get("A520", "2")
n, nz, dx = 48, 6, 0.05
L = n*dx
I, J = np.meshgrid(np.arange(n+1), np.arange(n+1), indexing="ij")
x, y = I*dx, J*dx
SK = float(os.environ.get("SKEW", "1"))
X = x + SK*0.1*y + SK*0.04*np.sin(2*np.pi*x/L)*np.sin(np.pi*y/L)
Y = y + SK*0.05*x + SK*0.03*np.sin(np.pi*x/L)*np.sin(2*np.pi*y/L)
Z = 0.06*np.exp(-((X-1.5)**2+(Y-0.9)**2)/0.05)
ROT = int(os.environ.get("ROT", "1"))
XCTRL = os.environ.get("XCTRL", "").replace(";", "\n")
control = "B 1 %g\nB 10 0.0 %g 0.0 %g 0.0 %g\n\nM 10 NP\n\nR 20 1\n" % (dx, L, L, nz*dx)
ctrl = ("A 10 5\nA 520 %s\nF 60 0.3\nF 72 0.9 1.5 0.5 1.1 0.34\nG 50 1\nG 52 1\nN 47 0.3\nM 10 NP\nP 10 0\nW 22 -9.81\n%s\n" % (a520, XCTRL))

def layout(nb):
    m = n//nb
    offs, blocks = [], []
    for bj in range(nb):
        for bi in range(nb):
            k = ((bi + 2*bj) % 4) if ROT else 0
            if bi == 0 and bj == 0:
                k = 0
            offs.append((bi*m, bj*m, k, m))
    for (i0, j0, k, m) in offs:
        blocks.append(tuple(np.rot90(A[i0:i0+m+1, j0:j0+m+1], k=-k).copy() for A in (X, Y, Z)))
    return offs, blocks

def gidx(i0, j0, k, m, a, b):
    q = np.rot90(np.arange(m*m).reshape(m, m), k=-k)[a, b]
    return i0 + q//m, j0 + q % m

ref = None
for nb in nbs:
    offs, blocks = layout(nb)
    npr = 1 if nb == 1 else np_
    d = os.path.join(S4, "runs/t3_nb%d_np%d" % (nb, npr))
    c = ctrl + ("G 1 1\n" if nb == 1 else "")
    rc = run4(d, control.replace("NP", str(npr)), c.replace("NP", str(npr)), blocks, np_=npr, steps=steps)
    if rc != 0:
        print(nb, "rc", rc, log_lines(d, "!!!")[:3])
    G = amr_grids(d, npr)
    if nb == 1:
        ref = level0(G, n, n)
        its = [int(l.split()[1]) for l in log_lines(d, "piter:")]
        print("1 block: |P| %.1f, pressure iterations per solve (REEFMG) %.1f" % (np.abs(ref["P"]).max(), np.mean(its) if its else -1))
        continue
    err = {f: 0.0 for f in FIELDS}
    m0 = n//nb
    for g in G:
        i0, j0, k, m = offs[g["blk"]]
        for f in FIELDS:
            for a in range(g["nx"]):
                for b in range(g["ny"]):
                    if g["id"] == 0:
                        if a >= m0 or b >= m0:
                            continue
                        I0, J0 = a+g["I0"], b+g["J0"]
                    else:
                        I0, J0 = gidx(i0, j0, k, m, a, b)
                    err[f] = max(err[f], np.abs(g[f][a, b]-ref[f][I0, J0]).max())
    its = [int(l.split()[4]) for l in log_lines(d, "NHFLOW AMR pressure: iterations")]
    print("%2d blocks: pressure iterations per solve %.1f (max %d); vs one block %s" % (nb*nb, np.mean(its), max(its), "  ".join("%s %.1e" % (f, err[f]) for f in FIELDS)))
