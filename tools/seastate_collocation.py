#!/usr/bin/env python3
# Architect: Hans Bihs
"""
REEF3D::SEASTATE online collocation (A 787 1): observations -> seastate-collocation.csv

Reads CSV files of observations with the columns time, lon, lat, hs (the altimeter.csv and insitu_*.csv of a
download kit, or any file with these columns) and writes the file a run reads with A 787 1:

    time,x,y,source,id,hs
    2022-01-29 06:12:31.000,-12345.6,234567.8,alti,jason-3,7.21

with the position in the model frame (the same frame options as seastate_forcing.py and met_forcing.py).

  source   --source NAME for all files, otherwise from the file name: altimeter* -> alti, insitu_* -> buoy,
           else the file name
  id       the column mission, platform, station or id if present, otherwise the part of the file name after
           'insitu_' (the station), otherwise empty
  filters  --start / --end (YYYY-MM-DD[ HH:MM:SS]), --domain X0 X1 Y0 Y1 (model frame), --hs-range LO HI,
           --coast D --bathy seastate-bathy.dat: drop points closer than D [m] to a dry node of the raster
           (z > -0.5 m, or outside it) for the sources of --coast-sources (default alti)

Example (validation 76, Malik):
  python3 seastate_collocation.py kit/storms/malik/altimeter.csv kit/storms/malik/insitu_*.csv \\
          --origin 3 57 --start "2022-01-27 12:00" --domain -400000 560000 -720000 740000 \\
          --coast 12000 --bathy seastate-bathy.dat -o seastate-collocation.csv
"""
import argparse
import csv
import datetime as dt
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from seastate_forcing import Frame   # noqa: E402

MISSIONS = ["cryosat-2", "jason-1", "jason-2", "jason-3", "saral", "sentinel-3a", "envisat", "topex", "ers-1", "ers-2",
            "sentinel-3b", "sentinel-6a"]


def parse_time(s):
    s = s.strip().replace("T", " ").rstrip("Z")
    for f in ("%Y-%m-%d %H:%M:%S.%f", "%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M", "%Y-%m-%d"):
        try:
            return dt.datetime.strptime(s, f)
        except ValueError:
            pass
    if len(s) >= 8 and s[:8].isdigit():
        v = float(s)
        d = int(v)
        hms = int(round((v - d) * 1e6))
        return dt.datetime(d // 10000, d // 100 % 100, d % 100, hms // 10000, hms // 100 % 100, hms % 100)
    raise ValueError("time '%s'" % s)


def read_bathy(path):
    with open(path) as f:
        lines = [ln.split("$")[0].strip() for ln in f]
    lines = [ln for ln in lines]
    k = 1
    while not lines[k]:
        k += 1
    nx, ny = (int(v) for v in lines[k].split()[:2])
    x0, y0, dx, dy = (float(v) for v in lines[k + 1].split()[:4])
    z = np.array(" ".join(lines[k + 2:]).split(), float)[:nx * ny].reshape(ny, nx)
    return x0, y0, dx, dy, z


def coast_distance(bathy, x, y):
    x0, y0, dx, dy, z = bathy
    wet = np.pad(z < -0.5, 1, constant_values=False)        # outside the raster: dry
    try:
        from scipy.ndimage import distance_transform_edt
        d = distance_transform_edt(wet, sampling=(dy, dx))
    except ImportError:
        jj, ii = np.nonzero(~wet)
        dry = np.column_stack([ii * dx, jj * dy])
        d = np.zeros(wet.shape)
        J, I = np.nonzero(wet)
        P = np.column_stack([I * dx, J * dy])
        for a in range(0, len(P), 2000):
            q = P[a:a + 2000]
            d[J[a:a + 2000], I[a:a + 2000]] = np.sqrt(((q[:, None, :] - dry[None, :, :]) ** 2).sum(-1)).min(1)
    i = np.clip(np.rint((x - x0) / dx).astype(int) + 1, 0, wet.shape[1] - 1)
    j = np.clip(np.rint((y - y0) / dy).astype(int) + 1, 0, wet.shape[0] - 1)
    return d[j, i]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+")
    ap.add_argument("--origin", type=float, nargs=2, metavar=("LON0", "LAT0"))
    ap.add_argument("--crs")
    ap.add_argument("--offset", type=float, nargs=2)
    ap.add_argument("--rotation", type=float, default=0.0)
    ap.add_argument("--source")
    ap.add_argument("--start")
    ap.add_argument("--end")
    ap.add_argument("--domain", type=float, nargs=4, metavar=("X0", "X1", "Y0", "Y1"))
    ap.add_argument("--hs-range", type=float, nargs=2, default=[0.05, 25.0])
    ap.add_argument("--coast", type=float, default=0.0)
    ap.add_argument("--bathy")
    ap.add_argument("--coast-sources", nargs="+", default=["alti"])
    ap.add_argument("-o", "--output", default="seastate-collocation.csv")
    a = ap.parse_args(argv)
    fr = Frame(a)
    t0 = parse_time(a.start) if a.start else None
    t1 = parse_time(a.end) if a.end else None
    bathy = read_bathy(a.bathy) if a.coast > 0.0 else None
    if a.coast > 0.0 and not a.bathy:
        sys.exit("--coast needs --bathy")

    rows = []
    for path in a.files:
        base = os.path.basename(path)
        src = a.source or ("alti" if base.startswith("altimeter") else "buoy" if base.startswith("insitu_") else os.path.splitext(base)[0])
        fid = os.path.splitext(base)[0][7:] if base.startswith("insitu_") else ""
        with open(path, newline="") as f:
            rd = csv.DictReader(f)
            cols = {c.lower(): c for c in rd.fieldnames}
            for need in ("time", "lon", "lat", "hs"):
                if need not in cols:
                    sys.exit("%s: no column %s" % (path, need))
            idcol = next((cols[c] for c in ("mission", "platform", "station", "id") if c in cols), None)
            T, LON, LAT, HS, ID = [], [], [], [], []
            for r in rd:
                try:
                    t = parse_time(r[cols["time"]])
                    lon, lat, hs = float(r[cols["lon"]]), float(r[cols["lat"]]), float(r[cols["hs"]])
                except (ValueError, TypeError):
                    continue
                if not (np.isfinite(hs) and a.hs_range[0] < hs < a.hs_range[1]):
                    continue
                if (t0 and t < t0) or (t1 and t > t1):
                    continue
                i = r[idcol].strip() if idcol else fid
                if i.isdigit() and src == "alti" and int(i) < len(MISSIONS):
                    i = MISSIONS[int(i)]
                T.append(t); LON.append(lon); LAT.append(lat); HS.append(hs); ID.append(i)
        if not T:
            continue
        x, y = fr.to_model(np.array(LON), np.array(LAT))
        keep = np.ones(len(T), bool)
        if a.domain:
            keep &= (x >= a.domain[0]) & (x < a.domain[1]) & (y >= a.domain[2]) & (y < a.domain[3])
        if bathy is not None and src in a.coast_sources:
            keep &= coast_distance(bathy, x, y) >= a.coast
        for k in np.flatnonzero(keep):
            rows.append((T[k], x[k], y[k], src, ID[k], HS[k]))
        print("%s: %d of %d observations (%s)" % (path, int(keep.sum()), len(T), src))

    rows.sort(key=lambda r: r[0])
    with open(a.output, "w") as f:
        f.write("time,x,y,source,id,hs\n")
        for t, x, y, s, i, hs in rows:
            f.write("%s,%.1f,%.1f,%s,%s,%.3f\n" % (t.strftime("%Y-%m-%d %H:%M:%S.%f")[:-3], x, y, s, i, hs))
    print("%s: %d observations" % (a.output, len(rows)))


if __name__ == "__main__":
    main()
