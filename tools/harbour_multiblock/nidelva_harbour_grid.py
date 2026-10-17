# P5d: two-block grid of the lower Nidelva and the harbour of Trondheim (G 50 1, DIVEMesh R 20 1) from a geo.dat
# raster (x y z per line, the local frame of georef.json; water z <= 0, land z > 0):
#   block 0 (the native grid): a curvilinear corridor along the river from the upstream end of the raster (inflow,
#            C 11 1 with B 60 1, W 10 Q) to the straight line y = Y0 where it enters the harbour block
#   block 1: the harbour and the fjord, a uniform box (DX) with i to the north and j to the west (right-handed); its
#            -i face (y = Y0) holds the river mouth of block 0 (point matching, the corridor ends on the box nodes),
#            its +i face (y = YN, in the fjord) is open: the water level (G 54 1 2 WL)
# The raster has no water depths: the bed is SCHEMATIC, from the distance to the shore: the river 4 m deep in
# the middle (D = 4 (1 - exp(-d/12 m))), harbour and fjord 8 m at 30 m from the quays and deeper offshore
# (8 (1 - exp(-d/8 m)) + 0.04 max(d - 60 m, 0), at most DMAX), smoothed over 15 m.  The vertical datum is moved
# up by WL (the still water level, default 10 m): z = WL - depth in the water, WL + min(z, 5 m) on land.
# usage: python3 nidelva_harbour_grid.py geo.dat outdir [--dx 10] [--q 100] [--wl 10] [--nz 5]
import argparse, os, json, numpy as np
from scipy import ndimage
from scipy.interpolate import RegularGridInterpolator, splprep, splev
from skimage.graph import route_through_array

ap = argparse.ArgumentParser()
ap.add_argument("geo")
ap.add_argument("out")
ap.add_argument("--dx", type=float, default=10.0, help="spacing of the harbour block and along the river (m)")
ap.add_argument("--nj", type=int, default=10, help="cells across the river corridor")
ap.add_argument("--box", type=float, nargs=4, default=[-400.0, 2300.0, -300.0, 2000.0], help="harbour block xmin xmax ymin ymax")
ap.add_argument("--mouth", type=float, nargs=2, default=[1180.0, 1280.0], help="river corridor on the south face of the box: x from to")
ap.add_argument("--q", type=float, default=100.0, help="river discharge (m3/s)")
ap.add_argument("--wl", type=float, default=10.0, help="still water level (m above the raster datum)")
ap.add_argument("--dmax", type=float, default=25.0, help="largest depth of the fjord (m)")
ap.add_argument("--nz", type=int, default=5, help="sigma layers (NHFLOW)")
ap.add_argument("--split", type=int, default=2, help="harbour blocks side by side along j (west-east; for the load of the ranks)")
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)

# ------------------------------------------------------------------ raster, water mask, schematic bed
d = np.loadtxt(a.geo)
xs, ys = np.unique(d[:, 0]), np.unique(d[:, 1])
h = xs[1]-xs[0]
Z = np.full((len(xs), len(ys)), np.nan)
Z[np.searchsorted(xs, d[:, 0]), np.searchsorted(ys, d[:, 1])] = d[:, 2]
Z = np.nan_to_num(Z, nan=10.0)
W = Z <= 0.0
lab, n = ndimage.label(~W)                       # land specks in the water (single raster points of the map)
size = ndimage.sum(~W, lab, range(1, n+1))
W |= np.isin(lab, 1+np.where(size < 40)[0])
lab, n = ndimage.label(W)                        # isolated ponds
size = ndimage.sum(W, lab, range(1, n+1))
W = np.isin(lab, 1+np.where(size >= 400)[0])
dist = ndimage.distance_transform_edt(W)*h       # distance to the shore (m)

X2, Y2 = np.meshgrid(xs, ys, indexing="ij")
river = W & (Y2 < a.box[2]+1.0)                  # the river south of the box (and the lower reach in the box below)
lower = W & (X2 > 1100.0) & (X2 < 1900.0) & (Y2 < 900.0) & (Y2 >= a.box[2])
dep_r = 4.0*(1.0-np.exp(-dist/12.0))
dep_h = np.minimum(8.0*(1.0-np.exp(-dist/8.0)) + 0.04*np.maximum(dist-60.0, 0.0), a.dmax)
dep = np.where(river | lower, dep_r, dep_h)*W
dep = ndimage.gaussian_filter(dep, 15.0/h)*W
bed = np.where(W, a.wl - np.maximum(dep, 0.05), a.wl + np.minimum(np.maximum(Z, 0.5), 5.0))
bed_at = RegularGridInterpolator((xs, ys), bed, bounds_error=False, fill_value=a.wl+5.0)

# ------------------------------------------------------------------ block 1: harbour and fjord
x0, x1, y0, y1 = a.box
DX = a.dx
nI, nJ = int(round((y1-y0)/DX)), int(round((x1-x0)/DX))
I, J = np.meshgrid(np.arange(nI+1), np.arange(nJ+1), indexing="ij")
BX, BY = x1 - DX*J, y0 + DX*I                    # i north, j west
BZ = bed_at(np.c_[BX.ravel(), BY.ravel()]).reshape(BX.shape)

# ------------------------------------------------------------------ block 0: river corridor
# centreline: least-cost path along the middle of the river from its upstream end to the mouth on y = Y0
m = river & (Y2 < y0)
dr = ndimage.distance_transform_edt(m)
cost = np.where(m, 1.0/(1.0+dr)**2, 1.0e6)
js = np.argmin(np.where(m.any(0), ys, 1e9))      # the southernmost water row: the upstream end
js += 3
cols = np.where(m[:, js])[0]
i_up = cols[np.argmin(abs(cols-np.mean(cols)))]
row = np.searchsorted(ys, y0-h)
cols = np.where(m[:, row])[0]
i_dn = cols[np.argmin(abs(xs[cols]-0.5*(a.mouth[0]+a.mouth[1])))]
path, _ = route_through_array(cost, (i_up, js), (i_dn, row), fully_connected=True, geometric=True)
P = np.array([(xs[p], ys[q]) for p, q in path])[::4]
tck, u = splprep(P.T, s=len(P)*8.0**2)
uu = np.linspace(0, 1, 4000)
cx, cy = splev(uu, tck)
s = np.r_[0, np.cumsum(np.hypot(np.diff(cx), np.diff(cy)))]
# the last LS m bend smoothly into a straight line to the north ending in the middle of the mouth
LS = 60.0
k = np.searchsorted(s, s[-1]-LS)
xm = 0.5*(a.mouth[0]+a.mouth[1])
cut = s[k]
L = cut + np.hypot(xm-cx[k], y0-cy[k])
ni = int(round(L/DX))
nj = a.nj
sn = np.linspace(0, L, ni+1)
c = np.c_[np.interp(sn, s, cx), np.interp(sn, s, cy)]
straight = np.c_[xm + 0.0*sn, y0 - (L-sn)]
wgt = np.clip((sn-(cut-LS))/(L-(cut-LS)), 0, 1)
wgt = wgt*wgt*(3-2*wgt)
c = (1-wgt)[:, None]*np.where((sn <= s[-1])[:, None], c, straight) + wgt[:, None]*straight
c[-1] = [xm, y0]
t = np.gradient(c, axis=0)
t /= np.hypot(t[:, 0], t[:, 1])[:, None]
nrm = np.c_[-t[:, 1], t[:, 0]]                   # left of the flow: j increases to the left (right-handed)
# curvature: the half-width on the inner side of a bend stays below 0.7 of the radius
dt_ = np.gradient(t, axis=0)/np.gradient(sn)[:, None]
kap = ndimage.gaussian_filter1d(t[:, 0]*dt_[:, 1] - t[:, 1]*dt_[:, 0], 2.0)   # > 0: turning left

# half-widths to the banks along the normals (water of the river), a dry margin of DX, smoothed along the river
wat = RegularGridInterpolator((xs, ys), m.astype(float), bounds_error=False, fill_value=0.0)
def extent(k, sign):
    r = np.arange(0.0, 250.0, 2.0)
    pts = c[k][None, :] + sign*r[:, None]*nrm[k][None, :]
    w = wat(pts) > 0.5
    out = np.where(~w)[0]
    return r[out[0]] if len(out) else r[-1]
hl = np.array([extent(k, 1.0) for k in range(ni+1)]) + DX
hr = np.array([extent(k, -1.0) for k in range(ni+1)]) + DX
hl = ndimage.gaussian_filter1d(np.maximum(hl, 2.0*DX), 4.0)
hr = ndimage.gaussian_filter1d(np.maximum(hr, 2.0*DX), 4.0)
R = 1.0/np.maximum(abs(kap), 1e-6)
hl = np.where(kap > 0, np.minimum(hl, 0.7*R), hl)
hr = np.where(kap < 0, np.minimum(hr, 0.7*R), hr)

# the last 15 cells blend into the straight cross-section of the mouth (nodes of the box)
NB = 15
wgt = np.clip((np.arange(ni+1)-(ni-NB))/NB, 0, 1)**2
hm = 0.5*(a.mouth[1]-a.mouth[0])
hl = (1-wgt)*hl + wgt*hm
hr = (1-wgt)*hr + wgt*hm
tt = np.linspace(-1.0, 1.0, nj+1)
CX = np.empty((ni+1, nj+1)); CY = np.empty((ni+1, nj+1))
for k in range(ni+1):
    off = np.where(tt < 0, tt*hr[k], tt*hl[k])
    CX[k] = c[k, 0] + off*nrm[k, 0]
    CY[k] = c[k, 1] + off*nrm[k, 1]
# the mouth: exactly the box nodes (j of the corridor to the west, as j of the box)
CX[-1] = a.mouth[1] - DX*np.arange(nj+1)
CY[-1] = y0
assert abs((a.mouth[1]-a.mouth[0]) - nj*DX) < 1e-9, "the mouth spans nj cells of the box"

# elliptic smoothing of the interior (Laplace on the nodes, boundaries fixed): no folded cells at the inner bends
for it in range(300):
    CX[1:-1, 1:-1] = 0.25*(CX[2:, 1:-1]+CX[:-2, 1:-1]+CX[1:-1, 2:]+CX[1:-1, :-2])
    CY[1:-1, 1:-1] = 0.25*(CY[2:, 1:-1]+CY[:-2, 1:-1]+CY[1:-1, 2:]+CY[1:-1, :-2])

def jac(X, Y):
    ax, ay = X[1:, :-1]-X[:-1, :-1], Y[1:, :-1]-Y[:-1, :-1]
    bx, by = X[:-1, 1:]-X[:-1, :-1], Y[:-1, 1:]-Y[:-1, :-1]
    return ax*by - ay*bx
for name, (X, Y) in (("river", (CX, CY)), ("harbour", (BX, BY))):
    jj = jac(X, Y)
    print("%-8s %4d x %3d cells, area %.1f .. %.1f m2, %s" % (name, X.shape[0]-1, X.shape[1]-1, jj.min(), jj.max(),
          "right-handed" if jj.min() > 0 else "FOLDED"))
CZ = bed_at(np.c_[CX.ravel(), CY.ravel()]).reshape(CX.shape)

# ------------------------------------------------------------------ output
def write_blocks(path, blocks):
    with open(path, "w") as f:
        f.write("%d\n" % len(blocks))
        for X, Y, Zb in blocks:
            f.write("%d %d\n" % (X.shape[0]-1, X.shape[1]-1))
            for j in range(X.shape[1]):
                for i in range(X.shape[0]):
                    f.write("%.6f %.6f %.4f\n" % (X[i, j], Y[i, j], Zb[i, j]))
# the harbour as SPLIT blocks along j (east to west), block 1 holds the mouth (j = 0 the east edge)
cuts = np.linspace(0, nJ, a.split+1).round().astype(int)
assert cuts[1] >= int(round((x1-a.mouth[0])/DX)), "the mouth lies in the first harbour block"
hb = [(BX[:, c0:c1+1], BY[:, c0:c1+1], BZ[:, c0:c1+1]) for c0, c1 in zip(cuts[:-1], cuts[1:])]
write_blocks(os.path.join(a.out, "blocks.dat"), [(CX, CY, CZ)] + hb)
np.savez_compressed(os.path.join(a.out, "grid.npz"), CX=CX, CY=CY, CZ=CZ, BX=BX, BY=BY, BZ=BZ, xs=xs, ys=ys, bed=bed, W=W)
H = a.wl + 5.0
open(os.path.join(a.out, "control.txt"), "w").write(
    "C 11 1\nC 12 21\nC 13 21\nC 14 21\nC 15 21\nC 16 3\n\nB 1 %g\nB 10 %g %g %g %g 0.0 %g\n\nM 10 NP\n\nR 20 1\n"
    % (H/a.nz, min(CX.min(), x0), max(CX.max(), x1), min(CY.min(), y0), y1, H))
json.dump({"river_cells": [ni, nj], "harbour_cells": [nI, nJ], "harbour_blocks": a.split, "river_length_m": L, "Q": a.q, "WL": a.wl,
           "mouth_j_of_box": [int(round((x1-a.mouth[1])/DX)), int(round((x1-a.mouth[0])/DX))],
           "bed": "schematic (distance to the shore)"}, open(os.path.join(a.out, "grid.json"), "w"), indent=1)
print("river corridor %.0f m, %d x %d cells; harbour %d x %d cells; mouth j %d..%d of the box" %
      (L, ni, nj, nI, nJ, round((x1-a.mouth[1])/DX), round((x1-a.mouth[0])/DX)))
