# P3 T9: mesh refinement with wetting and drying in the patches (G 30 1): a hump of water running up a
# plane beach (bed rising in x, the upper part dry), a static patch over the shoreline: Cartesian (bed
# by DIVEMesh S 61), curvilinear on the same nodes and turned by 90 deg (bed in curvgrid.dat): level 0
# per cell index and every patch, volume
import sys, glob, numpy as np
from p3lib import *
import p1lib
from readvts import readvts
steps = int(sys.argv[1]) if len(sys.argv) > 1 else 40
n, nz, dx = 64, 6, 0.0625
L = n*dx
control = "C 11 21\nC 12 21\nC 13 21\nC 14 21\nC 15 21\nC 16 3\n\nB 2 %d %d %d\nB 10 0.0 %g 0.0 %g 0.0 0.35\n\nS 61 2.0 4.0 0.0 4.0 0.0 0.3\n\nM 10 1\n" % (n, n, nz, L, L)
ctrl = "A 10 5\nA 520 %s\nF 60 0.2\nN 47 0.3\nN 48 0\nN 49 0.004\nM 10 1\nP 10 0\nP 30 0.04\nW 22 -9.81\nG 1 1\nG 4 8\nG 30 1\n" % os.environ.get("A520","2")
I, J = np.meshgrid(np.arange(n+1), np.arange(n+1), indexing="ij")
X0, Y0 = I*dx, J*dx
Zbed = np.clip((X0-2.0)/2.0*0.3, 0.0, None)          # S 61 2.0 4.0 ...: 0 at x = 2, 0.3 at x = 4
OX, OY = 1000.0, 2000.0
Xr, Yr = OX - (Y0-2.0), OY + (X0-2.0)
def boxes(rot):
    if not rot:
        return "F 72 1.0 1.5 1.5 2.5 0.3\nG 10 2.25 3.5 1.0 3.0\n"
    return "F 72 %.4f %.4f %.4f %.4f 0.3\nG 10 %.4f %.4f %.4f %.4f\n" % (OX-0.5, OX+0.5, OY-1.0, OY-0.5, OX-1.0, OX+1.0, OY+0.25, OY+1.5)
res = {}
for name, grid, rot in [("cart", None, False), ("curv_id", (X0, Y0, Zbed), False), ("rot90", (Xr, Yr, Zbed), True)]:
    c = ctrl + boxes(rot) + ("G 50 2\n" if grid is not None else "")
    d = os.path.join(S3, "runs/t9_"+name)
    rc = p1lib.run(d, control if grid is None else control.replace("S 61 2.0 4.0 0.0 4.0 0.0 0.3\n",""), c, grid=grid, steps=steps)
    st = final_state(d)
    log = [l.split() for l in open(os.path.join(d, "REEF3D_NHFLOW_AMR", "REEF3D_NHFLOW_AMR_log.dat")) if l.strip() and l.split()[0].isdigit()]
    last = sorted(glob.glob(os.path.join(d, "REEF3D_NHFLOW_AMR", "REEF3D-NHFLOW-AMR-*.vtm")))[-1]
    stamp = os.path.basename(last).split("-")[-1].split(".")[0]
    pf = sorted(f for f in glob.glob(os.path.join(d, "REEF3D_NHFLOW_AMR", "REEF3D-NHFLOW-AMR-%s-*" % stamp)) if not f.endswith(".vtm") and "-L0-" not in f)
    res[name] = (st, log, [readvts(f)[0] for f in pf])
    wet = st["WL"] > 0.002
    print("%-8s rc=%d patches %s  wet %d of %d  volume rel %s  max |U| %.3f" % (name, rc, log[-1][3], wet.sum(), wet.size, max(abs(float(l[8])) for l in log), np.abs(st["U"]).max()), p1lib.log_lines(d, "!!!"))
print("bed curv_id vs cart: %.1e" % np.abs(res["curv_id"][0]["bed"]-res["cart"][0]["bed"]).max())
for name, ref in [("curv_id","cart"), ("rot90","cart"), ("rot90","curv_id")]:
    st, log, pats = res[name]
    cst, clog, cpats = res[ref]
    U, V = st["U"], st["V"]
    cu, cv = (cst["U"], cst["V"]) if name=="curv_id" else (-cst["V"], cst["U"])
    pe = 0.0
    if len(pats)==len(cpats):
        for a, b in zip(pats, cpats):
            pe = max(pe, np.abs(a["eta"]-b["eta"]).max(), np.abs(a["u_top"]-(b["u_top"] if name=="curv_id" else -b["v_top"])).max())
    print("%-8s vs %-7s: eta %.1e  U %.1e  V %.1e  WL %.1e | patches %s %.1e" % (name, ref, np.abs(st["eta"]-cst["eta"]).max(), np.abs(U-cu).max(), np.abs(V-cv).max(), np.abs(st["WL"]-cst["WL"]).max(), "same" if len(pats)==len(cpats) else "COUNT", pe))
