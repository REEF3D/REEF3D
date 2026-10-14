# P5a T2: river reach over an erodible bed with a mound (bed on the nodes of curvgrid.dat): discharge inflow
# (B 60 1), stage outflow (F 62), bed friction; curvilinear grid on the Cartesian nodes, turned by 90 deg,
# turned by 30 deg at UTM coordinates, and a skewed (wavy) grid: bed per cell index against the grid on the
# nodes, bed volume (sum of bed change times cell area) against the sediment flux through the in- and outflow.
# environment: SOLVER nhflow|sflow, A520 (NHFLOW), XCTRL extra ctrl lines (;), NP ranks, MOUND height, SLOPE
# (steep mound for the sand slide: MOUND=0.08 SLOPE=0.15)
import sys, glob, numpy as np
from p5lib import *
steps = int(sys.argv[1]) if len(sys.argv) > 1 else 60
nx, ny, nz, dx = 40, 12, 6, 0.1
L, W, h = nx*dx, ny*dx, 0.25
np_ = int(os.environ.get("NP", "1"))
solver = os.environ.get("SOLVER", "nhflow")
bc = "C 11 1\nC 12 21\nC 13 21\nC 14 2\nC 15 21\nC 16 3\n"
ctl = control(nx, ny, nz, L, W, 0.5, bc=bc, np_=np_)
if solver == "nhflow":
    base = "A 10 5\nA 519 1\nA 520 %s\nB 50 0.01\n" % os.environ.get("A520", "0")
else:
    base = "A 10 2\nA 210 3\nA 211 4\nA 220 %s\nA 240 1\nB 50 0.01\n" % os.environ.get("A220", "0")
sed = "S 10 1\nS 11 2\nS 20 0.0005\nS 15 3\nS 17 1.0\nS 35 50.0\nS 41 1\nS 43 1\nS 44 1\nS 31 1\n"
base += "B 60 1\nW 10 %g\nF 60 %g\nF 62 %g\nN 47 0.3\nM 10 %d\nP 10 0\nW 22 -9.81\n" % (0.4*W*h, h, h, np_) + sed + os.environ.get("XCTRL", "").replace(";", "\n")
mound, slope = float(os.environ.get("MOUND", "0.04")), float(os.environ.get("SLOPE", "0.3"))
I, J = np.meshgrid(np.arange(nx+1), np.arange(ny+1), indexing="ij")
s0, n0 = I*dx, J*dx
Z = np.maximum(mound - np.hypot(s0-1.5, n0-0.5*W)*slope, 0.0)    # cone over the bed, off the centre line
X0, Y0 = 2.7e5, 7.04e6
res = {}
for name, rot in [("curv_id", "id"), ("rot90", "rot90"), ("rot30", "rot"), ("wavy", "wavy")]:
    if rot == "wavy":
        X, Y = nodes(nx, ny, dx, "id", wavy=0.25)
    elif rot == "rot":
        X, Y = nodes(nx, ny, dx, "rot", X0, Y0, 30.0)
    else:
        X, Y = nodes(nx, ny, dx, rot, 1000.0, 2000.0)
    d = os.path.join(S5, "runs/t2_%s_%s" % (solver, name))
    rc = run(d, ctl, base + "G 50 2\n", grid=(X, Y, Z), np_=np_, steps=steps)
    st = final_state(d)
    st0 = read_state(sorted(glob.glob(os.path.join(d, "dump", "state_*_r0.bin")))[0])
    # cell areas from the nodes
    A = 0.5*((X[1:,1:]-X[:-1,:-1])*(Y[:-1,1:]-Y[1:,:-1]) - (Y[1:,1:]-Y[:-1,:-1])*(X[:-1,1:]-X[1:,:-1]))
    dvol = ((st["bed"]-st0["bed"]).reshape(nx, ny)*A).sum()
    res[name] = st
    print("%-8s rc=%d t=%.3f  bed change min %.3e max %.3e  volume change %.3e  max|U| %.3f" % (name, rc, st["simtime"],
          (st["bed"]-st0["bed"]).min(), (st["bed"]-st0["bed"]).max(), dvol, np.abs(st.get("U", st.get("P"))).max()), log_lines(d, "!!!")[:2])
for name in ["rot90", "rot30"]:
    r, c = res[name], res["curv_id"]
    db = np.abs(r["bed"]-c["bed"]).max()
    print("%-8s vs curv_id: bed %.1e (rel to the change %.1e)  eta %.1e" % (name, db, db/max(np.abs(c["bed"]-Z.max()*0).max(),1e-30),
          np.abs(r["eta"]-c["eta"]).max()))
