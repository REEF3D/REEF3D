#!/usr/bin/env python3
# slide_terrain.py: the terrain and the initial slide mass of an SFLOW (or CFD) landslide case from two DEMs,
# one before and one after the event (GeoTIFF, same projected coordinate system, e.g. EPSG:25833).
#
#   python3 slide_terrain.py <pre.tif> <post.tif> <case dir> --bbox xmin ymin xmax ymax --dx 5
#                            --release xmin ymin xmax ymax [--thr 1.0] [--rotate deg]
#
# - the domain: --bbox (projected coordinates), resampled to the cell size --dx (block mean of the DEM cells)
# - the release area: inside the box --release, where the terrain fell by more than --thr metres
# - bed (geo.dat, for DIVEMesh G 10 1): the pre-event terrain, in the release area the post-event surface (the
#   sliding surface with the material left on it, as the reports measure it)
# - the initial level (sflow-level.dat, for REEF3D SFLOW F 73 1): the pre-event terrain in the release area
# - the reference (reference.npz): post - pre on the same grid (deposit > 0, depletion < 0), the release mask, the
#   origin of the local coordinates (x - xmin, y - ymin) for the evaluation
# All files in local coordinates (x - xmin, y - ymin, z unchanged).
import argparse, os, sys
import numpy as np


def read_dem(path, bbox, dx):
    import rasterio
    from rasterio.enums import Resampling
    from rasterio.windows import from_bounds
    xmin, ymin, xmax, ymax = bbox
    nx = int(round((xmax - xmin)/dx)); ny = int(round((ymax - ymin)/dx))
    with rasterio.open(path) as src:
        win = from_bounds(xmin, ymin, xmax, ymax, transform=src.transform)
        z = src.read(1, window=win, out_shape=(ny, nx), resampling=Resampling.average, masked=True, boundless=True,
                     fill_value=src.nodata if src.nodata is not None else -9999.0)
        crs = src.crs
    z = np.ma.filled(z.astype(float), np.nan)[::-1, :]     # row 0 = ymin
    z[z < -1000.0] = np.nan
    return z, nx, ny, crs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pre"); ap.add_argument("post"); ap.add_argument("case")
    ap.add_argument("--bbox", nargs=4, type=float, required=True)
    ap.add_argument("--release", nargs=4, type=float, required=True)
    ap.add_argument("--dx", type=float, default=5.0)
    ap.add_argument("--thr", type=float, default=1.0)
    a = ap.parse_args()
    xmin, ymin, xmax, ymax = a.bbox
    pre, nx, ny, crs = read_dem(a.pre, a.bbox, a.dx)
    post, _, _, _ = read_dem(a.post, a.bbox, a.dx)
    if np.isnan(pre).any() or np.isnan(post).any():
        print("warning: %d / %d cells without data (pre / post), filled with the other DEM"
              % (np.isnan(pre).sum(), np.isnan(post).sum()))
        pre = np.where(np.isnan(pre), post, pre); post = np.where(np.isnan(post), pre, post)
    xc = xmin + (np.arange(nx) + 0.5)*a.dx; yc = ymin + (np.arange(ny) + 0.5)*a.dx
    X, Y = np.meshgrid(xc, yc)
    rx0, ry0, rx1, ry1 = a.release
    inbox = (X >= rx0) & (X <= rx1) & (Y >= ry0) & (Y <= ry1)
    dz = post - pre
    release = inbox & (dz < -a.thr)
    bed = np.where(release, post, pre)
    level = np.where(release, pre, -9999.0)
    cell = a.dx*a.dx
    os.makedirs(a.case, exist_ok=True)
    with open(os.path.join(a.case, "geo.dat"), "w") as f:
        for j in range(ny):
            for i in range(nx):
                f.write("%.2f %.2f %.3f\n" % (xc[i] - xmin, yc[j] - ymin, bed[j, i]))
    with open(os.path.join(a.case, "sflow-level.dat"), "w") as f:
        f.write("%d %d %.3f %.3f %.3f %.3f\n" % (nx, ny, xc[0] - xmin, yc[0] - ymin, a.dx, a.dx))
        for j in range(ny):
            f.write(" ".join("%.3f" % v for v in level[j]) + "\n")
    np.savez_compressed(os.path.join(a.case, "reference.npz"), x=xc - xmin, y=yc - ymin, pre=pre, post=post, dz=dz,
                        release=release, origin=np.array([xmin, ymin]), dx=a.dx)
    dep = (~inbox) & (dz > a.thr)
    print("crs %s, grid %d x %d, dx %.1f m" % (crs, nx, ny, a.dx))
    print("release: %.3f km2, mobilised (pre - post) %.3f Mm3, max depth %.1f m"
          % (release.sum()*cell/1e6, -dz[release].sum()*cell/1e6, -dz[release].min()))
    print("outside the release box: deposit (post - pre > thr) %.3f Mm3 over %.3f km2, erosion %.3f Mm3"
          % (dz[dep].sum()*cell/1e6, dep.sum()*cell/1e6, -dz[(~inbox) & (dz < -a.thr)].sum()*cell/1e6))
    print("terrain %.1f - %.1f m; DIVEMesh: B 10 0 %.1f 0 %.1f <zmin> <zmax>, G 10 1" % (np.nanmin(bed), np.nanmax(pre), xmax - xmin, ymax - ymin))


if __name__ == "__main__":
    main()
