# P3 T7: NHFLOW with mesh refinement on curvilinear grids. Dam break (F 72 hump) in a closed basin with
# patches following the surface (G 20, G 21, regrid every step): Cartesian, curvilinear on the same
# nodes, turned by 90 deg: level 0 per cell index and every patch (AMR output) to round-off, mass.
# VAR: "" | sub (G 7 1) | vref (G 6 1) | static (G 10 box); NP ranks; A520
import sys, glob, numpy as np
from p3lib import *
import p1lib
from readvts import readvts
steps = int(sys.argv[1]) if len(sys.argv) > 1 else 30
NP = int(os.environ.get("NP","1"))
var = os.environ.get("VAR","")
n, nz, dx = 40, 8, 0.05
control = "C 11 21\nC 12 21\nC 13 21\nC 14 21\nC 15 21\nC 16 3\n\nB 2 %d %d %d\nB 10 0.0 %g 0.0 %g 0.0 0.5\n\nM 10 %d\n" % (n, n, nz, n*dx, n*dx, NP)
amr = {"": "G 1 1\nG 20 0.0005\nG 21 0.0001\nG 4 4\nG 2 1\nG 5 0\n",
       "sub": "G 1 1\nG 20 0.0005\nG 21 0.0001\nG 4 4\nG 2 1\nG 5 0\nG 7 1\n",
       "vref": "G 1 1\nG 20 0.0005\nG 21 0.0001\nG 4 4\nG 2 1\nG 5 0\nG 6 1\n",
       "two": "G 1 2\nG 20 0.0005\nG 21 0.0001\nG 4 4\nG 2 1\nG 5 0\n",
       "static": "G 1 1\nG 10 0.6 1.4 0.6 1.4\nG 4 4\n",
       "place": "G 1 1\nG 20 0.0005\nG 21 0.0001\nG 4 4\nG 2 1\nG 5 0\nG 40 1\n",
       "cdiff": "G 1 1\nG 20 0.0005\nG 21 0.0001\nG 4 4\nG 2 1\nG 5 0\nA 512 2\nW 2 0.0005\n"}[var]
base = os.environ.get("XCTRL","").replace(";","\n") + "A 10 5\nA 520 %s\nF 60 0.5\nN 47 0.5\nN 48 0\nN 49 0.004\nM 10 %d\nP 10 0\nP 30 0.04\nW 22 -9.81\n" % (os.environ.get("A520","2"), NP) + amr
I, J = np.meshgrid(np.arange(n+1), np.arange(n+1), indexing="ij")
X0, Y0 = I*dx, J*dx
Z = np.zeros_like(X0)
OX, OY = 1000.0, 2000.0
Xr, Yr = OX - (Y0-1.0), OY + (X0-1.0)            # turned by 90 deg about (1, 1), shifted
def box(rot):
    # the hump and two AMR gauges (P 51) at the same physical points of the turned frame
    if not rot:
        return "F 72 0.91 1.11 0.91 1.11 0.54\nP 51 1.013 1.013\nP 51 1.31 0.79\n"
    return "F 72 %.6f %.6f %.6f %.6f 0.54\nP 51 %.6f %.6f\nP 51 %.6f %.6f\n" % (OX-0.11, OX+0.09, OY-0.09, OY+0.11, OX-0.013, OY+0.013, OX+0.21, OY+0.31)
def gbox(rot):
    if var!="static":
        return ""
    return "" if not rot else ""
res = {}
for name, grid, rot in [("cart", None, False), ("curv_id", (X0, Y0, Z), False), ("rot90", (Xr, Yr, Z), True)]:
    c = base + box(rot) + ("G 50 2\n" if grid is not None else "")
    if var=="static" and rot:
        c = c.replace("G 10 0.6 1.4 0.6 1.4", "G 10 %.4f %.4f %.4f %.4f" % (OX-0.4, OX+0.4, OY-0.4, OY+0.4))
    d = os.path.join(S3, "runs/t7_%s_%s%d" % (var or "adapt", name, NP))
    rc = p1lib.run(d, control, c, grid=grid, np_=NP, steps=steps)
    st = [final_state(d, r) for r in range(NP)]
    log = [l.split() for l in open(os.path.join(d, "REEF3D_NHFLOW_AMR", "REEF3D_NHFLOW_AMR_log.dat")) if l.strip() and l.split()[0].isdigit()]
    last = sorted(glob.glob(os.path.join(d, "REEF3D_NHFLOW_AMR", "REEF3D-NHFLOW-AMR-*.vtm")))[-1]
    stamp = os.path.basename(last).split("-")[-1].split(".")[0]
    pf = sorted(f for f in glob.glob(os.path.join(d, "REEF3D_NHFLOW_AMR", "REEF3D-NHFLOW-AMR-%s-*" % stamp)) if not f.endswith(".vtm"))
    patches = [readvts(f)[0] for f in pf]
    gz = np.loadtxt(os.path.join(d, "REEF3D_NHFLOW_AMR", "REEF3D_NHFLOW_AMR_gauges.dat"))
    res[name] = (st, log, patches)
    res[name+"_g"] = gz
    print("%-8s rc=%d steps %s  patches %s cells %s  mass rel %s  mean piter %.1f  patch files %d" % (name, rc, log[-1][0], log[-1][3], log[-1][4], log[-1][8], np.mean([float(l[5]) for l in log]), len(patches)))
    for l in p1lib.log_lines(d, "!!!"):
        print("   ", l)
for name, ref in [("curv_id","cart"), ("rot90","cart"), ("rot90","curv_id")]:
    st, log, patches = res[name]
    cst, clog, cpat = res[ref]
    out = []
    for f in ["eta", "U", "V", "W", "P"]:
        e = 0.0
        for r in range(NP):
            a, b = st[r][f], cst[r][f]
            if name == "rot90" and f == "U": b = -cst[r]["V"]
            if name == "rot90" and f == "V": b = cst[r]["U"]
            e = max(e, np.abs(a-b).max())
        out.append("%s %.1e" % (f, e))
    pe = 0.0
    if len(patches) == len(cpat):
        for a, b in zip(patches, cpat):
            ue = a["u_top"] - (-b["v_top"] if name=="rot90" else b["u_top"])
            pe = max(pe, np.abs(a["eta"]-b["eta"]).max(), np.abs(ue).max())
    print("%-8s vs %-7s: level 0 %s | patches %s eta,u %.1e | layouts %s" % (name, ref, "  ".join(out), "same count" if len(patches)==len(cpat) else "COUNT %d vs %d" % (len(patches), len(cpat)), pe, "same" if [l[3:5] for l in log]==[l[3:5] for l in clog] else "DIFFER"))
g1, g2 = res["rot90_g"], res["curv_id_g"]
print("AMR gauges (P 51, the cell value on curvilinear grids) rot90 vs curv_id: eta %.1e, levels %s" % (np.abs(g1[:,1:3]-g2[:,1:3]).max(), "same" if (g1[:,3:]==g2[:,3:]).all() else "DIFFER"))
