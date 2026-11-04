#!/usr/bin/env python3
# Architect: Hans Bihs
"""
REEF3D::SFLOW maxima of a run (P 112, P 115): reading, points for ctrl.txt, comparison with observations.

  ctrl      a CSV of points (high-water marks, gauges) -> P 115 lines for ctrl.txt; lon/lat columns are
            projected with the frame of the run (W 205 lon0 lat0 rot, W 206 X0 Y0, as met_forcing.py), x/y columns
            are taken as they are
  read      REEF3D_SFLOW_MAX/REEF3D-SFLOW-MAX.pvtp -> .npz (cell centres x, y and the fields) or .csv
  info      the fields of a REEF3D-SFLOW-MAX.pvtp: range, highest water level and where
  compare   REEF3D-SFLOW-MAX-points.dat against observed peak levels (the same order as the P 115 lines):
            flooded points, bias, rms, r2; optional offset added to the model (e.g. a seasonal level the run leaves out)

The fields (per cell; NaN where a value does not exist):
  max_eta [m, relative to the still water level F 60], time_max_eta [s], max_depth [m], max_speed [m/s],
  time_max_speed [s], flood_duration [s], first_flood_time [s], inundation (0 never wet, 1 wet at P 113,
  2 dry at P 113 and flooded later), max_wind10 [m/s] and min_air_pressure [Pa] with W 200.

Examples:
  python3 sflow_max.py ctrl hwm.csv --lon lon --lat lat --frame -85.5 29.5 0 230000 190000 >> ctrl.txt
  python3 sflow_max.py read REEF3D_SFLOW_MAX/REEF3D-SFLOW-MAX.pvtp -o max.npz
  python3 sflow_max.py compare REEF3D_SFLOW_MAX/REEF3D-SFLOW-MAX-points.dat hwm.csv --obs peak_m --offset 0.27
"""
import argparse
import csv
import math
import os
import re
import sys

import numpy as np

FIELDS = ["max_eta", "time_max_eta", "max_depth", "max_speed", "time_max_speed", "flood_duration",
          "first_flood_time", "inundation", "max_wind10", "min_air_pressure"]


# ------------------------------------------------------------------ VTP (appended raw, as REEF3D writes it)
def read_vtp(path):
    raw = open(path, "rb").read()
    i = raw.index(b'<AppendedData encoding="raw">')
    j = raw.index(b"_", i) + 1
    head = raw[:i].decode("latin1")
    out = {}
    for typ, attrs, off in re.findall(r'<DataArray type="(\w+)"([^>]*)offset="(\d+)"', head):
        m = re.search(r'Name="([^"]+)"', attrs)
        name = m.group(1) if m else "points"
        nc = re.search(r'NumberOfComponents="(\d+)"', attrs)
        nc = int(nc.group(1)) if nc else 1
        off = int(off)
        n = int(np.frombuffer(raw[j + off:j + off + 4], dtype="<i4")[0])
        dt = {"Float32": "<f4", "Int32": "<i4", "Float64": "<f8"}[typ]
        a = np.frombuffer(raw[j + off + 4:j + off + 4 + n], dtype=dt)
        out[name] = a.reshape(-1, nc) if nc > 1 else a
    return out


def read_max(pvtp):
    """cell centres and the fields of all pieces of a REEF3D-SFLOW-MAX.pvtp"""
    folder = os.path.dirname(os.path.abspath(pvtp))
    pieces = re.findall(r'<Piece Source="([^"]+)"', open(pvtp).read())
    X, Y, F = [], [], {}
    for s in pieces:
        a = read_vtp(os.path.join(folder, s))
        con = a["connectivity"].reshape(-1, 3)[0::2]          # the first triangle of each cell
        P = a["points"]
        # nodes (i-1,j-1) and (i,j) of the first triangle: the diagonal of the cell
        X.append(0.5 * (P[con[:, 0], 0] + P[con[:, 2], 0]))
        Y.append(0.5 * (P[con[:, 0], 1] + P[con[:, 2], 1]))
        for k in FIELDS:
            if k in a:
                F.setdefault(k, []).append(a[k][0::2].astype(float))
    out = {"x": np.concatenate(X), "y": np.concatenate(Y)}
    out.update({k: np.concatenate(v) for k, v in F.items()})
    return out


# ------------------------------------------------------------------ commands
def frame(args):
    """lon/lat -> model coordinates: met_frame of REEF3D (met_forcing.GeoFrame, W 205 lon0 lat0 rot, W 206 X0 Y0)"""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from met_forcing import GeoFrame
    return GeoFrame(*args.frame).to_model


def cmd_ctrl(args):
    rows = list(csv.DictReader(open(args.csv)))
    if args.lon:
        if args.frame is None:
            sys.exit("lon/lat need --frame lon0 lat0 rot X0 Y0 (W 205, W 206 of the run)")
        f = frame(args)
    n = 0
    for r in rows:
        if args.lon:
            x, y = f(float(r[args.lon]), float(r[args.lat]))
        else:
            x, y = float(r[args.x]), float(r[args.y])
        print("P 115 %.2f %.2f" % (x, y))
        n += 1
    print("# %d points (sflow_max.py ctrl %s)" % (n, os.path.basename(args.csv)), file=sys.stderr)


def cmd_read(args):
    d = read_max(args.pvtp)
    if args.o.endswith(".npz"):
        np.savez_compressed(args.o, **d)
    else:
        keys = list(d)
        with open(args.o, "w") as o:
            o.write(",".join(keys) + "\n")
            for r in zip(*[d[k] for k in keys]):
                o.write(",".join("%.7g" % v for v in r) + "\n")
    print("%s: %d cells, fields %s" % (args.o, len(d["x"]), ", ".join(k for k in d if k not in ("x", "y"))))


def cmd_info(args):
    d = read_max(args.pvtp)
    print("%d cells" % len(d["x"]))
    for k in FIELDS:
        if k in d:
            v = d[k]
            ok = np.isfinite(v)
            print("  %-17s %8d values  %12.5g .. %12.5g" % (k, ok.sum(), np.nanmin(v) if ok.any() else np.nan,
                                                          np.nanmax(v) if ok.any() else np.nan))
    if "max_eta" in d and np.isfinite(d["max_eta"]).any():
        k = int(np.nanargmax(d["max_eta"]))
        print("highest water level %.3f m at x = %.1f, y = %.1f, t = %.0f s" % (d["max_eta"][k], d["x"][k], d["y"][k],
                                                                              d["time_max_eta"][k]))
    if "inundation" in d:
        print("cells flooded (dry at P 113, wet later): %d" % int(np.sum(d["inundation"] == 2)))


def read_points(path):
    names = None
    rows = []
    for l in open(path):
        if l.startswith("#"):
            continue
        if names is None:
            names = l.split()
            continue
        rows.append([float(v) for v in l.split()])
    a = np.array(rows)
    return {k: a[:, n] for n, k in enumerate(names)}


def cmd_compare(args):
    m = read_points(args.points)
    obs = np.array([float(r[args.obs]) if r[args.obs] not in ("", "nan") else np.nan
                    for r in csv.DictReader(open(args.csv))])
    if len(obs) != len(m["n"]):
        sys.exit("%d observations, %d points: the CSV must be the one of the P 115 lines" % (len(obs), len(m["n"])))
    mod = m["max_eta"] + args.offset
    ok = np.isfinite(mod) & np.isfinite(obs)
    e = mod[ok] - obs[ok]
    print("points %d, with an observation %d, wet in the model %d" % (len(obs), np.isfinite(obs).sum(), ok.sum()))
    if ok.sum() > 0:
        print("bias %+.3f m  rms %.3f m  max |error| %.3f m" % (e.mean(), np.sqrt(np.mean(e ** 2)), np.abs(e).max()))
    if ok.sum() > 2:
        print("r2 %.3f" % (np.corrcoef(mod[ok], obs[ok])[0, 1] ** 2))
    if args.o:
        with open(args.o, "w") as o:
            o.write("n,x,y,observed,model,error\n")
            for k in range(len(obs)):
                o.write("%d,%.2f,%.2f,%.4g,%.4g,%.4g\n" % (m["n"][k], m["x"][k], m["y"][k], obs[k], mod[k], mod[k] - obs[k]))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    a = sp.add_parser("ctrl")
    a.add_argument("csv")
    a.add_argument("--lon"); a.add_argument("--lat")
    a.add_argument("--x"); a.add_argument("--y")
    a.add_argument("--frame", type=float, nargs=5, metavar=("LON0", "LAT0", "ROT", "X0", "Y0"))
    a = sp.add_parser("read")
    a.add_argument("pvtp")
    a.add_argument("-o", default="max.npz")
    a = sp.add_parser("info")
    a.add_argument("pvtp")
    a = sp.add_parser("compare")
    a.add_argument("points")
    a.add_argument("csv")
    a.add_argument("--obs", required=True)
    a.add_argument("--offset", type=float, default=0.0)
    a.add_argument("-o")
    args = ap.parse_args()
    {"ctrl": cmd_ctrl, "read": cmd_read, "info": cmd_info, "compare": cmd_compare}[args.cmd](args)


if __name__ == "__main__":
    main()
