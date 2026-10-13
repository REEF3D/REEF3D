# P4 T5: mesh refinement inside and across the blocks: the four turned blocks of t2 with a static refinement
# box over the vertex (LEV levels) or patches following the surface (ADAPT=1: G 20, regrid every 4 steps):
# level 0 and every fine cell against one block with the same refinement (environment: A520, LEV, ADAPT,
# XCTRL extra ctrl lines separated by ;)
import sys, numpy as np
from p4lib import *

steps = int(sys.argv[1]) if len(sys.argv) > 1 else 20
np_ = int(sys.argv[2]) if len(sys.argv) > 2 else 1
schemes = os.environ.get("A520", "0 2").split()
LEV = int(os.environ.get("LEV", "1"))
ADAPT = int(os.environ.get("ADAPT", "0"))
n, nz, dx = 40, 6, 0.05
L = n*dx
h = n//2
I, J = np.meshgrid(np.arange(n+1), np.arange(n+1), indexing="ij")
x, y = I*dx, J*dx
SK = float(os.environ.get("SKEW", "1"))
X = x + SK*0.15*y + SK*0.04*np.sin(2*np.pi*x/L)*np.sin(np.pi*y/L)
Y = y + SK*0.05*x + SK*0.03*np.sin(np.pi*x/L)*np.sin(2*np.pi*y/L)
Z = 0.06*np.exp(-((X-1.3)**2+(Y-0.7)**2)/0.04)

def quad(i0, j0, k):
    return tuple(np.rot90(A[i0:i0+h+1, j0:j0+h+1], k=-k).copy() for A in (X, Y, Z))

def gidx(i0, j0, k, m, a, b):
    """global cell of cell (a,b) of the block quad(i0,j0,k) on a level with m x m cells per block"""
    q = np.rot90(np.arange(m*m).reshape(m, m), k=-k)[a, b]
    return i0 + q//m, j0 + q % m

offs = [(0, 0, 0), (h, 0, 1), (0, h, 2), (h, h, 3)]
blocks = [quad(*o) for o in offs]
amr = "G 1 %d\nG 3 2\nG 4 8\n" % LEV
amr += "G 2 4\nG 20 0.004\n" if ADAPT else "G 2 0\nG 10 0.55 1.75 0.6 1.6\n"
XCTRL = os.environ.get("XCTRL", "").replace(";", "\n")
control = "B 1 %g\nB 10 0.0 %g 0.0 %g 0.0 %g\n\nM 10 NP\n\nR 20 1\n" % (dx, L, L, nz*dx)

for a520 in schemes:
    ctrl = ("A 10 5\nA 509 0\nA 520 %s\nF 60 0.3\nF 72 0.65 1.45 0.65 1.45 0.34\nG 50 1\nG 52 1\nN 47 0.3\nM 10 NP\nP 10 0\nW 22 -9.81\n%s%s\n"
            % (a520, amr, XCTRL))
    R = {}
    for name, bl, npr in [("one", [(X, Y, Z)], 1), ("four", blocks, np_)]:
        d = os.path.join(S4, "runs/t5_%s_A%s_L%d_np%d" % (name, a520, LEV, npr))
        rc = run4(d, control.replace("NP", str(npr)), ctrl.replace("NP", str(npr)), bl, np_=npr, steps=steps)
        if rc != 0:
            print(name, "rc", rc, log_lines(d, "!!!")[:3])
        R[name] = amr_grids(d, npr)
    # the reference: fields per level and global fine cell
    ref = {}
    for G in R["one"]:
        l = G["lev"]
        for a in range(G["nx"]):
            for b in range(G["ny"]):
                ref[(l, G["I0"]+a, G["J0"]+b)] = {f: G[f][a, b] for f in FIELDS}
    v0 = {G["blk"]: (G["I0"], G["J0"]) for G in R["four"] if G["lev"] == 0 and G["id"] > 0}
    err = {}
    ncells = {}
    for G in R["four"]:
        l = G["lev"]
        s = 1 << l
        b = G["blk"]
        i0, j0, k = offs[b]
        for a in range(G["nx"]):
            for c in range(G["ny"]):
                Iv, Jv = G["I0"]+a, G["J0"]+c
                if b == 0:
                    if G["id"] == 0 and (Iv >= h or Jv >= h):
                        continue
                    gI, gJ = Iv, Jv
                else:
                    oa, ob = Iv - v0[b][0]*s, Jv - v0[b][1]*s
                    gI, gJ = gidx(i0*s, j0*s, k, h*s, oa, ob)
                key = (l, int(gI), int(gJ))
                if key not in ref:
                    err.setdefault(l, {}).setdefault("missing", 0)
                    err[l]["missing"] += 1
                    continue
                ncells[l] = ncells.get(l, 0) + 1
                for f in FIELDS:
                    e = np.abs(G[f][a, c]-ref[key][f]).max()
                    err.setdefault(l, {})[f] = max(err.get(l, {}).get(f, 0.0), e)
    nref = {}
    for key in ref:
        nref[key[0]] = nref.get(key[0], 0) + 1
    log = open(os.path.join(S4, "runs/t5_four_A%s_L%d_np%d" % (a520, LEV, np_), "REEF3D_NHFLOW_AMR", "REEF3D_NHFLOW_AMR_log.dat")).read().splitlines()
    m = [l.split() for l in log if l and l[0] != "#"]
    print("A520 %s, %d steps, %s: volume change %s" % (a520, steps, "adaptive" if ADAPT else "static box", m[-1][-2]))
    for l in sorted(err):
        print("  level %d: %d cells (one block %d)  %s" % (l, ncells.get(l, 0), nref.get(l, 0), "  ".join("%s %.1e" % (f, v) if f != "missing" else "missing %d" % v for f, v in err[l].items())))
