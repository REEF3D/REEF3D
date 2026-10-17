# P5d T1: open faces of the blocks 1.. (G 54): a river reach (skewed grid turned by 30 deg at (1000, 2000), a mound on
# the bed, Manning friction) with the discharge inflow at its -i end and the water level at its +i end, as one block
# (native in- and outflow, B 60 1, F 62) and split along i into blocks with the inflow and / or the outflow on a block
# 1.. (G 54 b 1 Q, G 54 b 2 level):
#   in:   block 1 the first NA cells (G 54 1 1 Q), block 0 the rest (native outflow)
#   out:  block 0 the first cells (native inflow), block 1 the last NA cells (G 54 1 2 h)
#   both: block 1 inflow, block 0 the middle, block 2 outflow (no in- or outflow on the native grid)
# per cell of the reach (I, J): eta, U, V (NHFLOW: the bottom layer) against one block
# environment: SOLVER sflow|nhflow, NP (ranks of the block runs), LAYOUTS, XCTRL (;), NORUN
import sys, os, re, glob, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../sediment_curvilinear"))
from p5lib import amr_grids2, log_lines, nodes, read_state
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../nhflow_multiblock"))
from p4lib import run4

T = os.path.dirname(os.path.abspath(__file__))
steps = int(sys.argv[1]) if len(sys.argv) > 1 else 200
solver = os.environ.get("SOLVER", "sflow")
np_ = int(os.environ.get("NP", "1"))
layouts = os.environ.get("LAYOUTS", "in out both").split()

nx, ny, nz, dx, NA = 48, 10, 4, 0.1, 16
L, W, h = nx*dx, ny*dx, 0.25
Q = 0.4*W*h
X, Y = nodes(nx, ny, dx, "rot", 1000.0, 2000.0, 30.0, wavy=0.25)
I, J = np.meshgrid(np.arange(nx+1), np.arange(ny+1), indexing="ij")
Z = np.maximum(0.08 - np.hypot(I*dx-2.4, J*dx-0.5*W)*0.15, 0.0)

if solver == "sflow":
    flow = "A 10 2\nA 210 3\nA 211 4\nA 220 0\nA 240 1\n"
else:
    flow = "A 10 5\nA 519 1\n"
ctrl0 = (flow + "B 50 0.01\nB 60 1\nW 10 %g\nF 50 2\nF 60 %g\nF 62 %g\nG 50 1\nN 47 0.3\nN 48 0\nN 49 0.005\nM 10 NP\nP 10 0\nW 22 -9.81\n"
         % (Q, h, h) + os.environ.get("XCTRL", "").replace(";", "\n") + "\n")

def sub(i0, i1):
    return tuple(A[i0:i1+1].copy() for A in (X, Y, Z))

# layout: (blocks as (first cell, last cell), C 11 / C 14 of the native grid, G 54 lines)
LAY = {
    "one":  ([(0, nx)], 1, 2, ""),
    "in":   ([(NA, nx), (0, NA)], 21, 2, "G 54 1 1 %g\n" % Q),
    "out":  ([(0, nx-NA), (nx-NA, nx)], 1, 21, "G 54 1 2 %g\n" % h),
    "both": ([(NA, nx-NA), (0, NA), (nx-NA, nx)], 21, 21, "G 54 1 1 %g\nG 54 2 2 %g\n" % (Q, h)),
}

def control(c11, c14, npr):
    return ("C 11 %d\nC 12 21\nC 13 21\nC 14 %d\nC 15 21\nC 16 3\n\nB 1 %g\nB 10 0.0 %g 0.0 %g 0.0 %g\n\nM 10 %d\n\nR 20 1\n"
            % (c11, c14, dx, L, W, nz*0.1, npr))

R = {}
LAY["one2"] = LAY["one"]
for name in ["one"] + layouts:
    bl, c11, c14, g54 = LAY[name]
    npr = 1 if name == "one" else np_
    d = os.path.join(T, "runs/t1_%s_%s_np%d" % (solver, name, npr))
    c = ctrl0.replace("NP", str(npr)) + g54 + ("G 1 1\n" if name.startswith("one") else "")
    rc = 0
    if not os.environ.get("NORUN"):
        rc = run4(d, control(c11, c14, npr), c, [sub(a, b) for a, b in bl], np_=npr, steps=steps)
    ring = [0, 0, 0, 0]
    for l in log_lines(d, "ring"):
        m = re.search(r"ring (\d+) (\d+) (\d+) (\d+)", l)
        if m:
            ring = [int(v) for v in m.groups()]
    t = {}
    if not os.path.exists(os.path.join(d, "amr", "amr_r0.bin")):
        # SFLOW on one block: no AMR dump, eta from the regression dump (the velocities only as face fluxes there)
        e = np.concatenate([np.asarray(read_state(sorted(glob.glob(os.path.join(d, "dump", "state_*_r%d.bin" % r)))[-1])["eta"])
                            for r in range(npr)]).reshape(nx, ny)
        for a in range(nx):
            for b in range(ny):
                t[(a, b)] = (e[a, b], np.nan, np.nan)
    for g in (amr_grids2(d, npr) if os.path.exists(os.path.join(d, "amr", "amr_r0.bin")) else []):
        if g["lev"] != 0:
            continue
        for a in range(g["nx"]):
            for b in range(g["ny"]):
                if g["blk"] == 0:
                    io, jo = g["I0"]+a-ring[0], g["J0"]+b-ring[2]
                    if io < 0 or io >= bl[0][1]-bl[0][0] or jo < 0 or jo >= ny:
                        continue
                    key = (bl[0][0]+io, jo)
                else:
                    key = None
                    # block b: its own cells from the first cell of its level-0 grid (I0 = vi)
                    key = (bl[g["blk"]][0]+a, b)
                t[key] = (g["eta"][a, b], g["U"][a, b], g["V"][a, b])
    R[name] = t
    print("%-5s rc=%d cells %d" % (name, rc, len(t)), log_lines(d, "!!!")[:2], log_lines(d, "Inflow block")[:1])

ref = R["one"]
print("one block: eta %.4f .. %.4f" % (min(v[0] for v in ref.values()), max(v[0] for v in ref.values())))
for name in layouts:
    t = R[name]
    com = [k for k in ref if k in t]
    e = [np.nanmax([abs(t[k][m]-ref[k][m]) for k in com]) if com and not np.isnan(ref[com[0]][m]) else float("nan") for m in range(3)]
    worst = max(com, key=lambda k: abs(t[k][0]-ref[k][0])) if com else None
    print("%-5s vs one (%d of %d cells): eta %.1e U %.1e V %.1e  (largest eta difference at %s)" % (name, len(com), len(ref), *e, worst))
if os.environ.get("PAIRS"):
    for a_ in layouts:
        for b_ in layouts:
            if a_ < b_:
                com = [k for k in R[a_] if k in R[b_]]
                print("%-5s vs %-5s: eta %.1e U %.1e" % (a_, b_, max(abs(R[a_][k][0]-R[b_][k][0]) for k in com),
                                                       np.nanmax([abs(R[a_][k][1]-R[b_][k][1]) for k in com])))
    E = np.zeros((nx, ny))
    for k in ref:
        if k in R[layouts[0]]:
            E[k] = R[layouts[0]][k][0]-ref[k][0]
    np.set_printoptions(precision=1, linewidth=200)
    print((E*1e6).T.round(0)[:, 14:34])
