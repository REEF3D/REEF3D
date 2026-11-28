#!/usr/bin/env python3
# valdez_terrain.py: the terrain, the initial slides and the observations of NTHMP landslide benchmark 7 (Port Valdez,
# Alaska, 27.3.1964) for the slide layer (W 150, W 161 1) under SFLOW or NHFLOW.
#
#   python3 valdez_terrain.py <benchmark_7_data dir> <case dir> --dx 40 [--slides hpv sbm] [--box x0 x1 y0 y1]
#                             [--zclip 150]
#
# The benchmark data (www1.udel.edu/kirby/landslide/problems/benchmark_7.html, benchmark_7_data.zip): ASCII grids
# 570 x 2316 on a regular longitude-latitude grid (about 13.4 x 16.5 m), WGS84
#   valdez.bathymetry.long / .latg   longitude, latitude of the nodes
#   valdez.1964.bathymetry           post-earthquake terrain (m, positive up, MLLW): the surface the slides move on
#   valdez.HPV_1964, valdez.SBM_1964 thickness of the head-of-bay slide (HPV64) and the Shoup Bay moraine slide (SBM64)
#   *.shp                            observed inundation at the head of the bay, debris line of the first wave,
#                                    McKinley Street (lines, longitude-latitude)
# Local coordinates: x = R cos(61.1 deg) (lon - lon_min), y = R (lat - lat_min), metres (equirectangular; the scale
# error over the 9 km in latitude is 0.13 %).
# Writes: geo.dat (DIVEMesh G 10 1: the terrain, clipped at --zclip), slide-level.dat (W 161 1: terrain + thickness of
# the chosen slides, -9999 elsewhere), reference.npz (terrain, thicknesses, observations, lines; for rheo_eval.py
# valdez).
import argparse, os, struct
import numpy as np

R_EARTH = 6371000.0
LAT_REF = 61.1

# observations (benchmark 7 page; Plafker et al. 1969, Wilson & Torum 1972). Points given by the benchmark: the hotel
# and the navigation light. The others are placed by the description and the shoreline of the data (approximate).
OBS = [
    # name, lon, lat, quantity, observed value, note
    ("Valdez Hotel", -146.26701, 61.115823, "depth", 0.46, "0.46 m in the hotel on McKinley St (second wave)"),
    ("Narrows light", -146.652749, 61.081169, "eta", 11.0, "pedestal 11 m above MLLW overtopped"),
]
# shore stretches for the highest runup (lon_min, lon_max, lat_min, lat_max): the largest water level on land in them
SHORES = [
    ("Anderson Bay (S shore)", (-146.600, -146.560, 61.076, 61.092), 24.0, "camp destroyed, runup 24 m"),
    ("Cliff Mine (N shore)", (-146.645, -146.605, 61.105, 61.130), 52.0, "driftwood 52 m, splash 67 m"),
    ("Jackson Point, Dayville (S shore)", (-146.460, -146.400, 61.080, 61.098), 9.5, "cannery, waves to 9.5 m"),
    ("Valdez waterfront", (-146.300, -146.255, 61.105, 61.125), None, "inundation beyond McKinley St (about 300 m)"),
]


def read_shp(path):
    """polylines / polygons of an ESRI shapefile: a list of (n, 2) arrays (lon, lat)"""
    b = open(path, "rb").read()
    flen = struct.unpack(">i", b[24:28])[0]*2
    pos, parts_out = 100, []
    while pos < flen:
        clen = struct.unpack(">i", b[pos + 4:pos + 8])[0]*2
        rec = b[pos + 8:pos + 8 + clen]
        st = struct.unpack("<i", rec[:4])[0]
        if st in (3, 5, 13, 15, 23, 25):
            npart, npt = struct.unpack("<ii", rec[36:44])
            parts = list(struct.unpack("<%di" % npart, rec[44:44 + 4*npart])) + [npt]
            pts = np.frombuffer(rec[44 + 4*npart:44 + 4*npart + 16*npt], dtype="<f8").reshape(npt, 2)
            for k in range(npart):
                parts_out.append(pts[parts[k]:parts[k + 1]].copy())
        pos += 8 + clen
    return parts_out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("data"); ap.add_argument("case")
    ap.add_argument("--dx", type=float, default=40.0)
    ap.add_argument("--slides", nargs="+", default=["hpv", "sbm"], choices=["hpv", "sbm"])
    ap.add_argument("--box", nargs=4, type=float, default=None, help="x0 x1 y0 y1 (local metres) of the domain")
    ap.add_argument("--zclip", type=float, default=150.0)
    ap.add_argument("--hmin", type=float, default=0.05, help="slide thickness below which a cell has no slide")
    a = ap.parse_args()
    from scipy.interpolate import RegularGridInterpolator
    rd = lambda f: np.loadtxt(os.path.join(a.data, f))
    lon, lat = rd("valdez.bathymetry.long"), rd("valdez.bathymetry.latg")
    z = rd("valdez.1964.bathymetry")
    th = {"hpv": rd("valdez.HPV_1964"), "sbm": rd("valdez.SBM_1964")}
    lo, la = lon[0, :], lat[:, 0]
    lon0, lat0 = lo.min(), la.min()
    cx = np.deg2rad(1.0)*R_EARTH*np.cos(np.deg2rad(LAT_REF)); cy = np.deg2rad(1.0)*R_EARTH
    xs, ys = (lo - lon0)*cx, (la - lat0)*cy
    tolocal = lambda lon_, lat_: ((np.asarray(lon_) - lon0)*cx, (np.asarray(lat_) - lat0)*cy)
    X0, X1, Y0, Y1 = a.box if a.box else (0.0, xs.max(), 0.0, ys.max())
    nx, ny = int((X1 - X0)//a.dx), int((Y1 - Y0)//a.dx)
    xc = X0 + (np.arange(nx) + 0.5)*a.dx; yc = Y0 + (np.arange(ny) + 0.5)*a.dx
    YY, XX = np.meshgrid(yc, xc, indexing="ij")
    # block mean: the mean of the data nodes in each cell (the data are 13-16 m, the cells larger); bilinear below 25 m
    def resample(f):
        if a.dx >= 25.0:
            ix = np.floor((xs - X0)/a.dx).astype(int); iy = np.floor((ys - Y0)/a.dx).astype(int)
            IX, IY = np.meshgrid(ix, iy)
            ok = (IX >= 0) & (IX < nx) & (IY >= 0) & (IY < ny)
            s = np.zeros((ny, nx)); c = np.zeros((ny, nx))
            np.add.at(s, (IY[ok], IX[ok]), f[ok]); np.add.at(c, (IY[ok], IX[ok]), 1.0)
            if (c == 0).any():
                g = RegularGridInterpolator((ys, xs), f, bounds_error=False, fill_value=None)((YY, XX))
                return np.where(c > 0, s/np.maximum(c, 1), g)
            return s/c
        return RegularGridInterpolator((ys, xs), f, bounds_error=False, fill_value=None)((YY, XX))
    bed = resample(z)
    hs = {k: np.maximum(resample(v), 0.0) for k, v in th.items()}
    h = sum(hs[k] for k in a.slides)
    h[h < a.hmin] = 0.0
    bedc = np.minimum(bed, a.zclip)
    os.makedirs(a.case, exist_ok=True)
    xl, yl = xc - X0, yc - Y0
    with open(os.path.join(a.case, "geo.dat"), "w") as f:
        for j in range(ny):
            for i in range(nx):
                f.write("%.2f %.2f %.3f\n" % (xl[i], yl[j], bedc[j, i]))
    level = np.where(h > 0, bedc + h, -9999.0)
    with open(os.path.join(a.case, "slide-level.dat"), "w") as f:
        f.write("%d %d %.3f %.3f %.3f %.3f\n" % (nx, ny, xl[0], yl[0], a.dx, a.dx))
        for j in range(ny):
            f.write(" ".join("%.3f" % v for v in level[j]) + "\n")
    lines = {}
    for n in ("InundationLine", "InundationLineIslands", "InundationLineFirstWave", "McKinley_St"):
        p = os.path.join(a.data, n + ".shp")
        if os.path.exists(p):
            lines[n] = [np.column_stack(tolocal(q[:, 0], q[:, 1])) - [X0, Y0] for q in read_shp(p)]
    obs = [(n, *(np.array(tolocal(lo_, la_)) - [X0, Y0]), q, v) for n, lo_, la_, q, v, _ in OBS]
    shores = []
    for n, (l0, l1, b0, b1), v, _ in SHORES:
        (x0_, y0_), (x1_, y1_) = np.array(tolocal(l0, b0)) - [X0, Y0], np.array(tolocal(l1, b1)) - [X0, Y0]
        shores.append((n, x0_, x1_, y0_, y1_, np.nan if v is None else v))
    np.savez_compressed(os.path.join(a.case, "reference.npz"), x=xl, y=yl, dx=a.dx, bed=bed, h_hpv=hs["hpv"],
                        h_sbm=hs["sbm"], slides=np.array(a.slides), origin_lonlat=np.array([lon0, lat0]),
                        origin_xy=np.array([X0, Y0]), scale=np.array([cx, cy]),
                        obs=np.array([(o[1], o[2], o[4]) for o in obs]), obs_names=np.array([o[0] for o in obs]),
                        obs_qty=np.array([o[3] for o in obs]),
                        shores=np.array([s[1:] for s in shores]), shore_names=np.array([s[0] for s in shores]),
                        **{"line_%s_%d" % (k, i): q for k, v in lines.items() for i, q in enumerate(v)})
    cell = a.dx*a.dx
    print("grid %d x %d, dx %.1f m, %.1f x %.1f km" % (nx, ny, a.dx, nx*a.dx/1e3, ny*a.dx/1e3))
    for k in ("hpv", "sbm"):
        print("%s: %.1f Mm3 over %.2f km2%s" % (k.upper(), hs[k].sum()*cell/1e6, (hs[k] > a.hmin).sum()*cell/1e6,
              " (released)" if k in a.slides else ""))
    print("terrain %.1f - %.1f m (clipped at %.0f); DIVEMesh: B 2 %d %d <nz>, B 10 0 %.1f 0 %.1f %.0f %.0f, G 10 1"
          % (bed.min(), bed.max(), a.zclip, nx, ny, nx*a.dx, ny*a.dx, np.floor(bed.min()/10)*10 - 10, a.zclip + 10))
    for o in obs:
        j, i = int(o[2]//a.dx), int(o[1]//a.dx)
        if 0 <= i < nx and 0 <= j < ny:
            print("  %-15s x %8.1f y %7.1f terrain %6.1f m" % (o[0], o[1], o[2], bed[j, i]))


if __name__ == "__main__":
    main()
