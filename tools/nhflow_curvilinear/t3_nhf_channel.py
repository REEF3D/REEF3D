# P3 T3 (gate): a straight wave channel, linear waves kh ~ 2 generated at x- (Dirichlet, B 98 3) and
# absorbed at x+ (active beach, B 99 3): Cartesian, curvilinear turned by 30 deg at UTM coordinates,
# and on an orthogonal mapped grid (stretched along the channel, turned by 30 deg). Surface along the
# channel against the Cartesian run, Poisson iterations.
import sys, numpy as np, re
from p3lib import *
import p1lib
T_END = float(os.environ.get("TEND","8.0"))
nx, ny, nz, dx, h = 200, 4, 10, 0.1, 0.5
L, W = nx*dx, ny*dx
control = "C 11 6\nC 12 3\nC 13 3\nC 14 7\nC 15 21\nC 16 3\n\nB 2 %d %d %d\nB 10 0.0 %g 0.0 %g 0.0 %g\n\nM 10 1\nM 20 2\n" % (nx, ny, nz, L, W, h)
A520 = int(os.environ.get("A520","2"))
WG = os.environ.get("WG","dirichlet")
wg = "B 98 3\nB 99 3\n" if WG=="dirichlet" else "B 98 2\nB 99 1\n"
AMR = os.environ.get("AMR","")
amr = {"": "", "adapt": "G 1 1\nG 20 0.002\nG 4 4\nG 2 4\nP 30 1.0\n", "sub": "G 1 1\nG 20 0.002\nG 4 4\nG 2 4\nG 7 1\nP 30 1.0\n"}[AMR]
base = amr + ("A 10 5\nA 520 %d\nB 90 1\nB 92 1\nB 93 0.02 1.0\n" + wg + "F 60 %g\nN 41 %g\nN 47 0.5\nN 48 0\nN 49 0.01\nM 10 1\nP 10 0\nW 22 -9.81\n") % (A520, h, T_END)
X0, Y0 = 2.7e5, 7.04e6
th = np.radians(30.0)
I, J = np.meshgrid(np.arange(nx+1), np.arange(ny+1), indexing="ij")
xi = I/nx
s_unif = xi*L
s_map = L*(xi + 0.05*np.sin(2*np.pi*xi))       # stretched, monotone (min spacing 0.69 dx)
n_ = J*dx
def turned(s, n):
    return X0 + s*np.cos(th) - n*np.sin(th), Y0 + s*np.sin(th) + n*np.cos(th)
def zones(rot):
    # relaxation zones (WG relax): generation 3 m at the start, beach 5 m at the end, lines across the channel
    if WG=="dirichlet":
        return ""
    def P(sv, nv):
        if not rot:
            return sv, nv
        return X0 + sv*np.cos(th) - nv*np.sin(th), Y0 + sv*np.sin(th) + nv*np.cos(th)
    a, b = P(0.0, -1.0), P(0.0, W+1.0)
    c, e = P(L, -1.0), P(L, W+1.0)
    return "B 108 %.4f %.4f %.4f %.4f 3.0\nB 107 %.4f %.4f %.4f %.4f 5.0\n" % (a[0], b[0], a[1], b[1], c[0], e[0], c[1], e[1])
res = {}
for name, grid, sgrid in [("cart", None, s_unif), ("rot30", turned(s_unif, n_), s_unif), ("map0", (s_map, n_), s_map), ("map30", turned(s_map, n_), s_map)]:
    ctrl = base + zones(name in ("rot30","map30")) + ("G 50 2\nB 105 30.0 %.3f %.3f\n" % (X0, Y0) if name in ("rot30","map30") else ("G 50 2\n" if grid is not None else ""))
    d = os.path.join(S3, "runs/t3_"+WG+AMR+"_"+name)
    g = None if grid is None else (grid[0], grid[1], np.full(grid[0].shape, -h) if False else np.zeros(grid[0].shape))
    rc = p1lib.run(d, control, ctrl, grid=g)
    st = final_state(d)
    it = [int(m.group(1)) for m in (re.search(r"piter:\s*(\d+)", l) for l in open(os.path.join(d,"reef.log"), errors="replace")) if m]
    if AMR:
        lg = [l.split() for l in open(os.path.join(d, "REEF3D_NHFLOW_AMR", "REEF3D_NHFLOW_AMR_log.dat")) if l.strip() and l.split()[0].isdigit()]
        it = [float(l[5]) for l in lg]
        print("   %s: patches at the end %s, columns %s, volume rel %s" % (name, lg[-1][3], lg[-1][4], lg[-1][8]))
    eta = st["eta"].reshape(nx, ny)
    sc = 0.5*(sgrid[1:,0] + sgrid[:-1,0])          # cell centres along the channel
    res[name] = (sc, eta.mean(axis=1), np.mean(it) if it else 0.0, st["simtime"], rc)
    print("%-6s rc=%d t=%.3f mean piter %.2f  max |eta| %.4f" % (name, rc, st["simtime"], res[name][2], np.abs(eta).max()))
sc0, e0 = res["cart"][0], res["cart"][1]
for name, ref in [("rot30","cart"), ("map0","cart"), ("map30","cart"), ("map30","map0")]:
    sc0, e0 = res[ref][0], res[ref][1]
    sc, e = res[name][0], res[name][1]
    m = (sc > 1.0) & (sc < 14.0)
    ei = np.interp(sc[m], sc0, e0)
    print("%-6s vs %-5s: max |d eta| / H %.4f  rms / H %.4f  iterations %+.1f %%" % (name, ref, np.abs(e[m]-ei).max()/0.02, np.sqrt(np.mean((e[m]-ei)**2))/0.02, 100.0*(res[name][2]/res[ref][2]-1.0)))
np.save(os.path.join(S3, "runs/t3_profiles.npy"), np.array([res[k][0] for k in res] + [res[k][1] for k in res], dtype=object), allow_pickle=True)
