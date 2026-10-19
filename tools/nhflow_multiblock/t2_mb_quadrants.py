# P4 T2: four blocks around a regular vertex, their indices turned by 0, 90, 180, 270 deg, on a wavy skewed
# grid with a bumpy bed; NHFLOW dam break over the vertex: one block vs. four, per cell, all fields
# (environment: A520 the pressure schemes, SKEW 0 a Cartesian grid, XCTRL extra ctrl lines separated by ;)
import sys, numpy as np
from p4lib import *

steps = int(sys.argv[1]) if len(sys.argv) > 1 else 30
np_ = int(sys.argv[2]) if len(sys.argv) > 2 else 1
schemes = os.environ.get("A520", "0 1 2").split()
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
    """the quadrant (i0..i0+h, j0..j0+h) with its index frame turned by k quarter turns"""
    return tuple(np.rot90(A[i0:i0+h+1, j0:j0+h+1], k=-k).copy() for A in (X, Y, Z))

def gidx(i0, j0, k, a, b):
    """the global cell of cell (a,b) of the block quad(i0,j0,k)"""
    q = np.rot90(np.arange(h*h).reshape(h, h), k=-k)[a, b]
    return i0 + q//h, j0 + q % h

offs = [(0, 0, 0), (h, 0, 1), (0, h, 2), (h, h, 3)]
blocks = [quad(*o) for o in offs]
XCTRL = os.environ.get("XCTRL", "").replace(";", "\n")
control = "B 1 %g\nB 10 0.0 %g 0.0 %g 0.0 %g\n\nM 10 NP\n\nR 20 1\n" % (dx, L, L, nz*dx)

for a520 in schemes:
    ctrl = ("A 10 5\nA 520 %s\nF 60 0.3\nF 72 0.65 1.45 0.65 1.45 0.34\nG 50 1\nG 52 1\nN 47 0.3\nM 10 NP\nP 10 0\nW 22 -9.81\n%s\n"
            % (a520, XCTRL))
    R = {}
    for name, bl, npr in [("one", [(X, Y, Z)], 1), ("four", blocks, np_)]:
        d = os.path.join(S4, "runs/t2_%s_A%s_np%d" % (name, a520, npr))
        c = ctrl + ("G 1 1\n" if name == "one" else "")
        rc = run4(d, control.replace("NP", str(npr)), c.replace("NP", str(npr)), bl, np_=npr, steps=steps)
        if rc != 0:
            print(name, "rc", rc, log_lines(d, "!!!")[:3])
        R[name] = amr_grids(d, npr)
    ref = level0(R["one"], n, n)
    flds = FIELDS + tuple(f for f in ("EV", "KIN") if f in ref)
    err = {f: 0.0 for f in flds}
    for G in R["four"]:
        if G["blk"] == 0 and G["id"] == 0:
            for f in flds:
                i0, j0 = G["I0"], G["J0"]
                for a in range(G["nx"]):
                    for b in range(G["ny"]):
                        I0, J0 = i0+a, j0+b
                        if I0 < h and J0 < h:
                            err[f] = max(err[f], np.abs(G[f][a, b]-ref[f][I0, J0]).max())
            continue
        i0, j0, k = offs[G["blk"]]
        for f in flds:
            e = max(np.abs(G[f][a, b]-ref[f][gidx(i0, j0, k, a, b)]).max() for a in range(G["nx"]) for b in range(G["ny"]))
            err[f] = max(err[f], e)
    log = open(os.path.join(S4, "runs/t2_four_A%s_np%d" % (a520, np_), "REEF3D_NHFLOW_AMR", "REEF3D_NHFLOW_AMR_log.dat")).read().splitlines()
    m = [l.split() for l in log if l and l[0] != "#"]
    print("A520 %s, %d steps, t %.3f, |eta| %.3f |U| %.3f |P| %.1f: four blocks vs one: %s; volume change %s"
          % (a520, steps, R["one"][0]["simtime"], np.abs(ref["eta"]).max(), np.abs(ref["U"]).max(), np.abs(ref["P"]).max(),
             "  ".join("%s %.1e" % (f, err[f]) for f in flds), m[-1][-2]))
