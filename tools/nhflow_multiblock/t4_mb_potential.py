# P4 T4: potential flow initialisation (I 11, the var 44 solve) across blocks: a channel (block 0, discharge
# inflow at x = 0, stage outflow at x = L) blocked by a dry obstacle; the water goes around it through the
# region above the channel, split into NB blocks along x (their indices turned).  One block vs. 1 + NB
# blocks: the initial velocities, the iterations of the composite solve against the block count
# (environment: NBS the counts, STEPS time steps after the initialisation)
import sys, numpy as np
from p4lib import *

np_ = int(sys.argv[1]) if len(sys.argv) > 1 else 1
nbs = [int(v) for v in os.environ.get("NBS", "1 3 8 15").split()]
steps = int(os.environ.get("STEPS", "0"))
nx, ny0, ny1, nz, dx = 120, 12, 24, 6, 0.05
L, H = nx*dx, (ny0+ny1)*dx
I, J = np.meshgrid(np.arange(nx+1), np.arange(ny0+ny1+1), indexing="ij")
x, y = I*dx, J*dx
SK = float(os.environ.get("SKEW", "1"))
X = x + SK*0.03*np.sin(np.pi*y/H)*np.sin(2*np.pi*x/L)
Y = y + SK*0.04*np.sin(2*np.pi*x/L)*np.sin(np.pi*y/H)
xc = x
# dry: the obstacle in the channel (x 2.4 .. 3.6), the region above the channel left and right of the bypass
Z = np.where(((xc > 2.4) & (xc < 3.6) & (y < ny0*dx+1e-9)) | ((y > ny0*dx+1e-9) & ((xc < 1.2) | (xc > 4.8))), 0.6, 0.0)
control = "C 11 1\nC 12 21\nC 13 21\nC 14 2\nC 15 21\nC 16 3\n\nB 1 %g\nB 10 0.0 %g 0.0 %g 0.0 %g\n\nM 10 NP\n\nR 20 1\n" % (dx, L, H, nz*dx)
ctrl = ("A 10 5\nA 520 2\nB 60 1\nF 60 0.3\nF 62 0.3\nG 50 1\nG 52 1\nI 11 1\nN 47 0.3\nM 10 NP\nP 10 0\nW 10 0.05\nW 22 -9.81\n%s\n"
        % os.environ.get("XCTRL", "").replace(";", "\n"))

def layout(nb):
    if nb == 0:
        return [(0, 0, 0, nx, ny0+ny1)], [(X, Y, Z)]
    offs = [(0, 0, 0, nx, ny0)]
    w = nx//nb
    for b in range(nb):
        offs.append((b*w, ny0, (b % 3)+1, w, ny1))
    blocks = []
    for (i0, j0, k, mi, mj) in offs:
        blocks.append(tuple(np.rot90(A[i0:i0+mi+1, j0:j0+mj+1], k=-k).copy() for A in (X, Y, Z)))
    return offs, blocks

def gidx(i0, j0, k, mi, mj, a, b):
    q = np.rot90(np.arange(mi*mj).reshape(mi, mj), k=-k)[a, b]
    return i0 + q//mj, j0 + q % mj

ref = None
for nb in [0] + nbs:
    offs, blocks = layout(nb)
    npr = 1 if nb == 0 else np_
    d = os.path.join(S4, "runs/t4_nb%d_np%d" % (nb, npr))
    c = ctrl + ("G 1 1\n" if nb == 0 else "")
    rc = run4(d, control.replace("NP", str(npr)), c.replace("NP", str(npr)), blocks, np_=npr, steps=steps)
    if rc != 0:
        print(nb, "rc", rc, log_lines(d, "!!!|rror")[:3])
    G = amr_grids(d, npr)
    if nb == 0:
        ref = level0(G, nx, ny0+ny1)
        it = log_lines(d, "lapiter")
        print("one block: lapiter %s, max |U| %.4f |V| %.4f, U in the bypass %.4f" % (it[-1].split()[-1] if it else "?", np.abs(ref["U"]).max(), np.abs(ref["V"]).max(), np.abs(ref["U"][60, 24]).max()))
        continue
    err = {f: 0.0 for f in ("U", "V", "W")}
    for g in G:
        i0, j0, k, mi, mj = offs[g["blk"]]
        for f in err:
            for a in range(g["nx"]):
                for b in range(g["ny"]):
                    if g["id"] == 0:
                        if b >= ny0:
                            continue
                        I0, J0 = a+g["I0"], b+g["J0"]
                    else:
                        I0, J0 = gidx(i0, j0, k, mi, mj, a, b)
                    err[f] = max(err[f], np.abs(g[f][a, b]-ref[f][I0, J0]).max())
    it = log_lines(d, "potential flow on the blocks")
    print("%2d blocks: %s; vs one block %s" % (nb+1, it[-1].split(":")[-1] if it else "?", "  ".join("%s %.1e" % (f, err[f]) for f in err)))
