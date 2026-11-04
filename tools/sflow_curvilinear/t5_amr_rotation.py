# T5: mesh refinement on curvilinear grids: dam break with patches following the surface jump
# (G 20, regrid every 4 steps): Cartesian = curvilinear on the same nodes = grid rotated by 90 deg,
# per level-0 cell, for the hydrostatic, non-hydrostatic (A 220 2) and Boussinesq (A 220 4, a small
# step: undular bore) equations, with and without subcycling (G 7 1), on one and on two ranks with placement (G 40 1)
import sys, re, numpy as np
from p1lib import *

steps = int(sys.argv[1]) if len(sys.argv) > 1 else 40
n, dx = 120, 0.25
base = "A 10 2\nA 210 3\nA 211 4\nA 240 1\n%s\nF 60 1.0\nF 72 12.0 18.0 12.0 18.0 2.0\n\nN 41 100.0\nN 47 0.3\n\nM 10 %d\n\nW 22 -9.81\n\nG 1 2\nG 2 4\nG 3 2\nG 20 0.05\nG 4 8\n"
variants = [("hyd", "A 220 0\n", 1), ("hyd_sub", "A 220 0\nG 7 1\n", 1), ("nh", "A 220 2\n", 1),
            ("nh_sub", "A 220 2\nG 7 1\n", 1), ("bous", "A 220 4\n", 1), ("bous_sub", "A 220 4\nG 7 1\n", 1), ("bous_mpi2_place", "A 220 4\nG 40 1\n", 2), ("nh_mpi2_place", "A 220 2\nG 40 1\n", 2),
            ("hyd_box", "A 220 0\n", 1)]
if len(sys.argv) > 2:
    variants = [v for v in variants if v[0] in sys.argv[2:]]

I, J = np.meshgrid(np.arange(n+1), np.arange(n+1), indexing="ij")
X0, Y0 = I*dx, J*dx
Z = np.zeros_like(X0)
Xr, Yr = 1000.0-(Y0-15.0), 2000.0+(X0-15.0)

def merged(dd, np_):
    st = [final_state(dd, r) for r in range(np_)]
    return {k: np.concatenate([s[k] for s in st]) for k in ("eta", "WL")}, st[0]["simtime"]

def amrlog(dd):
    rows = [l.split() for l in open(os.path.join(dd, "REEF3D_SFLOW_AMR", "REEF3D_SFLOW_AMR_log.dat")) if l.strip() and l.split()[0].isdigit()]
    return rows

for name, key, np_ in variants:
    control = "B 1 %g\nB 10 0.0 %g 0.0 %g 0.0 0.25\n\nM 10 %d\nM 20 2\n" % (dx, n*dx, n*dx, np_)
    c = base % (key, np_)
    if name.startswith("bous"):
        # undular bore of a small step (the strong dam break makes the Boussinesq solution rough);
        # Krylov tolerance at round-off (the default N 44 gives agreement to about 1e-11)
        c = c.replace("F 72 12.0 18.0 12.0 18.0 2.0", "F 72 12.0 18.0 12.0 18.0 1.1").replace("G 20 0.05", "G 20 0.01") + "N 44 1e-14\n"
    if name == "hyd_box":
        c = c.replace("G 2 4\n", "G 2 0\n").replace("G 20 0.05\n", "G 10 8.0 22.0 8.0 22.0\n")
    R = {}
    for v, grid in [("cart", None), ("id", (X0, Y0, Z)), ("rot90", (Xr, Yr, Z))]:
        cc = c if v == "cart" else c + "G 50 2\n"
        if v == "rot90":
            cc = re.sub(r"F 72 12.0 18.0 12.0 18.0 (\S+)", r"F 72 997.0 1003.0 1997.0 2003.0 \1", cc)
        dd = os.path.join(S, "runs/t5_%s_%s" % (name, v))
        rc = run(dd, control, cc, grid=grid, np_=np_, steps=steps)
        if rc != 0:
            print(name, v, "FAILED rc", rc, log_lines(dd, "!!!")[:3]); R[v] = None; continue
        R[v] = merged(dd, np_) + (amrlog(dd),)
    if any(R[v] is None for v in R):
        continue
    lc, li, lr = R["cart"][2][-1], R["id"][2][-1], R["rot90"][2][-1]
    print("%-14s patches/cells %s/%s %s/%s %s/%s | eta curv=cart %.2e rot90 %.2e | mass rel. change cart %s id %s rot90 %s"
          % (name, lc[3], lc[4], li[3], li[4], lr[3], lr[4],
             np.abs(R["id"][0]["eta"]-R["cart"][0]["eta"]).max(), np.abs(R["rot90"][0]["eta"]-R["cart"][0]["eta"]).max(),
             lc[6], li[6], lr[6]))
