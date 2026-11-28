#!/usr/bin/env python3
# Architect: Hans Bihs
"""Terrain and an initial slide mound for the two-layer SFLOW model (W 150, W 161) from one DEM.

  slide_mound.py dem.tif <case dir> --bbox xmin ymin xmax ymax --dx 25
                 --mound x0 y0 L W theta_deg p V [--zmax 1200]

  dem.tif   the terrain (e.g. a post-event topo-bathy DEM with the scar), averaged onto the cells of size dx of the box
  --mound   a rotated super-Gaussian slide mass on the terrain, as the tsunami source of Shugar et al. (2026, Tracy Arm):
            h(x, y) = A exp(-(|x'/L|^p + |y'/W|^p)/2), x' = (x - x0) cos(theta) + (y - y0) sin(theta),
            y' = -(x - x0) sin(theta) + (y - y0) cos(theta), the amplitude A from the volume V [m^3]
  --zmax    the terrain is clipped at zmax (the mountains above the water and the slide do not matter)

Writes into <case dir>:
  geo.dat          x y z of the cell centres for DiveMESH (G 10 1, G 15 2), local coordinates x - xmin, y - ymin
  slide-level.dat  the top of the initial slide for SFLOW W 161 1 (the format of sflow-level.dat, F 73:
                   "nx ny x0 y0 dx dy", then the levels, x fastest, -9999 where no slide)
  reference.npz    pre (the terrain), dz (- the mound), release (the mound above 0.5 m), origin, dx: for
                   rheo_eval.py tsunamirunup
Needs numpy and rasterio.
"""
import os, argparse, math
import numpy as np


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dem"); ap.add_argument("case")
    ap.add_argument("--bbox", nargs=4, type=float, required=True)
    ap.add_argument("--dx", type=float, default=25.0)
    ap.add_argument("--mound", nargs=7, type=float, required=True, metavar=("x0", "y0", "L", "W", "theta", "p", "V"))
    ap.add_argument("--zmax", type=float, default=1.0e9)
    a = ap.parse_args()

    import rasterio
    from rasterio.enums import Resampling
    from rasterio.windows import from_bounds
    xmin, ymin, xmax, ymax = a.bbox
    nx, ny = int(round((xmax - xmin)/a.dx)), int(round((ymax - ymin)/a.dx))
    with rasterio.open(a.dem) as src:
        win = from_bounds(xmin, ymin, xmax, ymax, transform=src.transform)
        z = src.read(1, window=win, out_shape=(ny, nx), resampling=Resampling.average, masked=True, boundless=True)
        crs = src.crs
    z = np.ma.filled(z.astype(float), np.nan)[::-1, :]        # row 0 = ymin
    if np.isnan(z).any():
        print("warning: %d cells without data, set to the nearest value along the rows" % np.isnan(z).sum())
        for r in z:
            m = np.isnan(r)
            if m.any() and not m.all():
                r[m] = np.interp(np.where(m)[0], np.where(~m)[0], r[~m])
    z = np.minimum(z, a.zmax)

    x0, y0, L, W, th, p, V = a.mound
    th = math.radians(th)
    xc = xmin + (np.arange(nx) + 0.5)*a.dx; yc = ymin + (np.arange(ny) + 0.5)*a.dx
    X, Y = np.meshgrid(xc, yc)
    xp = (X - x0)*math.cos(th) + (Y - y0)*math.sin(th)
    yp = -(X - x0)*math.sin(th) + (Y - y0)*math.cos(th)
    I2 = 4.0*L*W*2.0**(2.0/p)/p**2*math.gamma(1.0/p)**2
    A = V/I2
    h = A*np.exp(-(np.abs(xp/L)**p + np.abs(yp/W)**p)/2.0)
    rel = h > 0.5
    level = np.where(rel, z + h, -9999.0)

    os.makedirs(a.case, exist_ok=True)
    with open(os.path.join(a.case, "geo.dat"), "w") as f:
        for j in range(ny):
            for i in range(nx):
                f.write("%.2f %.2f %.3f\n" % (xc[i] - xmin, yc[j] - ymin, z[j, i]))
    with open(os.path.join(a.case, "slide-level.dat"), "w") as f:
        f.write("%d %d %.3f %.3f %.3f %.3f\n" % (nx, ny, xc[0] - xmin, yc[0] - ymin, a.dx, a.dx))
        for j in range(ny):
            f.write(" ".join("%.3f" % v if v > -9000 else "-9999" for v in level[j]) + "\n")
    np.savez_compressed(os.path.join(a.case, "reference.npz"), x=xc - xmin, y=yc - ymin, pre=z, post=z, dz=-h,
                        release=rel, origin=np.array([xmin, ymin]), dx=a.dx)
    cell = a.dx*a.dx
    print("crs %s, grid %d x %d, dx %.1f m, terrain %.1f - %.1f m" % (crs, nx, ny, a.dx, z.min(), z.max()))
    print("mound: amplitude %.1f m, volume on the grid %.3f Mm3 (V %.3f Mm3), %d cells, top of the slide %.0f - %.0f m"
          % (A, h[rel].sum()*cell/1e6, V/1e6, rel.sum(), (z + h)[rel].min(), (z + h)[rel].max()))
    print("DIVEMesh: B 10 0 %.1f 0 %.1f %.0f %.0f, G 10 1" % (xmax - xmin, ymax - ymin, math.floor(z.min() - 10),
                                                         math.ceil(np.max(np.where(rel, z + h, z)) + 10)))


if __name__ == "__main__":
    main()
