#!/usr/bin/env python3
# Architect: Hans Bihs
"""
REEF3D atmospheric forcing (W 200 2): writer, reader and converters of the binary met file met-forcing.r3m.

  grid      10 m wind u10, v10 and mean sea level pressure on a lon/lat (or projected) grid from NetCDF, e.g. ERA5
            (u10, v10, msl), NORA3 or BARRA-C2 (uas, vas, psl) -> met-forcing.r3m on a regular grid in model
            coordinates; wind vectors turned from (east, north) to the model axes incl. the meridian convergence
  uniform   spatially uniform time series (nx = ny = 1) from a table: time U10 dir [pa]
  analytic  test fields: u10, v10, pa linear in x, y and t (exact under bilinear and linear interpolation)
  info      header and records of a met file

Model coordinates: as tools/seastate_forcing.py, --origin LON0 LAT0 (transverse Mercator on a sphere, the frame of
W 205 LON0 LAT0 ROT) or --crs EPSG:xxxx --offset X0 Y0 (pyproj), and --rotation ROT of the model x axis
counter-clockwise from east [deg]. The frame is stored in the file header and checked against W 205 / W 206.

File format version 1 (little-endian, src/met_file.h)
  header 144 bytes: char[8] "REEF3DMT"; int32 version, nx, ny, nvar, moving, nrec, frame, epsg;
                    float64 x0, y0, dx, dy, lon0, lat0, rot, X0, Y0, avg, reserved[3]
  records:          float64 time [s since 1970 UTC], xr0, yr0, reserved; float32 u10[ny][nx], v10[ny][nx], (pa[ny][nx])

Examples
  python3 met_forcing.py grid era5_2018_michael.nc --origin -85.5 29.5 --domain 0 400000 -200000 200000 \\
          --spacing 5000 -o met-forcing.r3m
  python3 met_forcing.py uniform series.txt -o met-forcing.r3m      # lines: YYYYMMDD.HHMMSS U10 dir [pa]
  python3 met_forcing.py analytic --domain 0 10000 0 5000 --spacing 500 --start 20180101.000000 --hours 6 --step 3600
  python3 met_forcing.py info met-forcing.r3m
"""
import argparse
import datetime as dt
import math
import os
import struct
import sys

import numpy as np

MAGIC = b"REEF3DMT"
VERSION = 1
HEADER = struct.Struct("<8s8i13d")          # 8 + 32 + 104 = 144 bytes
RECHEAD = struct.Struct("<4d")              # 32 bytes
EPOCH = dt.datetime(1970, 1, 1)


def datetime_to_epoch(t):
    return (t - EPOCH).total_seconds()


def stamp_to_epoch(v):
    """YYYYMMDD.HHMMSS (float or string) -> seconds since 1970-01-01 UTC, as met_datetime."""
    s = "%.6f" % float(v)
    d, f = s.split(".")
    t = dt.datetime(int(d[0:4]), int(d[4:6]), int(d[6:8])) + dt.timedelta(hours=int(f[0:2]), minutes=int(f[2:4]), seconds=int(f[4:6]))
    return datetime_to_epoch(t)


def epoch_to_str(t):
    return (EPOCH + dt.timedelta(seconds=float(t))).strftime("%Y-%m-%d %H:%M:%S")


def write_r3m(path, x0, y0, dx, dy, times, u, v, pa=None, frame=None, avg=0.0, origins=None):
    """times: seconds since 1970 (nrec); u, v (and pa): arrays (nrec, ny, nx) [m/s along the model axes, Pa];
    frame: dict(kind=1, lon0, lat0, rot, X0, Y0) or dict(kind=2, epsg, rot, X0, Y0) or None;
    origins: (nrec, 2) record origins for a moving grid, None: fixed grid."""
    u = np.asarray(u, dtype="<f4")
    v = np.asarray(v, dtype="<f4")
    nrec, ny, nx = u.shape
    nvar = 2 if pa is None else 3
    if pa is not None:
        pa = np.asarray(pa, dtype="<f4")
    times = np.asarray(times, float)
    if nrec != len(times) or v.shape != u.shape or (pa is not None and pa.shape != u.shape):
        raise ValueError("u, v, pa must be (nrec, ny, nx) with nrec = len(times)")
    if nrec > 1 and np.any(np.diff(times) <= 0.0):
        raise ValueError("the record times must increase")
    fr = frame or {}
    kind = int(fr.get("kind", 0))
    header = HEADER.pack(MAGIC, VERSION, nx, ny, nvar, 0 if origins is None else 1, nrec, kind, int(fr.get("epsg", 0)),
                         x0, y0, dx, dy, fr.get("lon0", 0.0), fr.get("lat0", 0.0), fr.get("rot", 0.0),
                         fr.get("X0", 0.0), fr.get("Y0", 0.0), avg, 0.0, 0.0, 0.0)
    with open(path, "wb") as out:
        out.write(header)
        for n in range(nrec):
            xr, yr = (x0, y0) if origins is None else origins[n]
            out.write(RECHEAD.pack(times[n], xr, yr, 0.0))
            out.write(u[n].tobytes())
            out.write(v[n].tobytes())
            if pa is not None:
                out.write(pa[n].tobytes())


def read_r3m(path):
    with open(path, "rb") as f:
        h = HEADER.unpack(f.read(HEADER.size))
        if h[0] != MAGIC:
            raise ValueError("%s: not a REEF3D met file" % path)
        version, nx, ny, nvar, moving, nrec, frame, epsg = h[1:9]
        x0, y0, dx, dy, lon0, lat0, rot, X0, Y0, avg = h[9:19]
        m = dict(version=version, nx=nx, ny=ny, nvar=nvar, moving=moving, nrec=nrec, frame=frame, epsg=epsg,
                 x0=x0, y0=y0, dx=dx, dy=dy, lon0=lon0, lat0=lat0, rot=rot, X0=X0, Y0=Y0, avg=avg)
        t, org, fields = [], [], []
        for n in range(nrec):
            r = RECHEAD.unpack(f.read(RECHEAD.size))
            t.append(r[0])
            org.append((r[1], r[2]))
            a = np.frombuffer(f.read(4 * nvar * nx * ny), dtype="<f4").reshape(nvar, ny, nx)
            fields.append(a)
        m["time"] = np.array(t)
        m["origin"] = np.array(org)
        a = np.array(fields)
        m["u"], m["v"] = a[:, 0], a[:, 1]
        m["pa"] = a[:, 2] if nvar == 3 else None
    return m


# ------------------------------------------------------------------ commands
def frame_of(args):
    """Frame object of tools/seastate_forcing.py and the header entries."""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import seastate_forcing as sf
    fr = sf.Frame(args)
    if args.crs:
        epsg = int(args.crs.split(":")[-1]) if ":" in args.crs else 0
        off = args.offset or [0.0, 0.0]
        hdr = dict(kind=2, epsg=epsg, rot=args.rotation, X0=off[0], Y0=off[1])
    else:
        hdr = dict(kind=1, lon0=args.origin[0], lat0=args.origin[1], rot=args.rotation, X0=0.0, Y0=0.0)
    return sf, fr, hdr


def first_var(nc, names):
    for n in names.split(","):
        if nc.has(n):
            return n
    return None


def grid(args):
    sf, fr, hdr = frame_of(args)
    nc = sf.NC(args.input)
    vu, vv = first_var(nc, args.u), first_var(nc, args.v)
    vp = first_var(nc, args.p) if args.p else None
    if vu is None or vv is None:
        sys.exit("no wind variables %s / %s in %s" % (args.u, args.v, args.input))
    u = nc.get(vu)
    v = nc.get(vv)
    pr = nc.get(vp) if vp else None
    lon = nc.get(args.lon)
    lat = nc.get(args.lat)
    times = sf.read_times(nc, args.time)
    keep = [k for k, t in enumerate(times) if (args.start is None or sf.stamp(t) >= args.start) and (args.end is None or sf.stamp(t) <= args.end)]
    if not keep:
        sys.exit("no times in the window")

    x0, x1, y0, y1 = args.domain
    nx = int(round((x1 - x0) / args.spacing)) + 1
    ny = int(round((y1 - y0) / args.spacing)) + 1
    X, Y = np.meshgrid(x0 + args.spacing * np.arange(nx), y0 + args.spacing * np.arange(ny))
    LON, LAT = fr.to_geo(X, Y)
    north = fr.north(LON, LAT)

    if lon.ndim == 1 and lat.ndim == 1:
        from scipy.interpolate import RegularGridInterpolator
        flip = lat[0] > lat[-1]
        la = lat[::-1] if flip else lat
        lonq = np.where(LON < lon.min(), LON + 360.0, LON) if lon.max() > 180.0 else LON

        def interp(field):
            fld = field[::-1, :] if flip else field
            g = RegularGridInterpolator((la, lon), fld, bounds_error=False, fill_value=None)
            return g(np.column_stack([LAT.reshape(-1), lonq.reshape(-1)])).reshape(LAT.shape)
    else:
        from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
        from scipy.spatial import Delaunay
        sx, sy = fr.to_model(lon.reshape(-1), lat.reshape(-1))
        pts = np.column_stack([sx, sy])
        tri = Delaunay(pts)

        def interp(field):
            vals = field.reshape(-1)
            a = LinearNDInterpolator(tri, vals)(X, Y)
            bad = np.isnan(a)
            if bad.any():
                a[bad] = NearestNDInterpolator(pts, vals)(X[bad], Y[bad])
            return a

    U = np.empty((len(keep), ny, nx), np.float32)
    V = np.empty_like(U)
    P = np.empty_like(U) if pr is not None else None
    for n, k in enumerate(keep):
        ue = interp(u[k])
        vn = interp(v[k])
        # (east, north) -> model frame: north at angle 'north', east at north - 90 deg
        U[n] = ue * np.sin(north) + vn * np.cos(north)
        V[n] = -ue * np.cos(north) + vn * np.sin(north)
        if P is not None:
            P[n] = interp(pr[k]) * args.p_scale
    T = [datetime_to_epoch(times[k]) for k in keep]
    write_r3m(args.output, x0, y0, args.spacing, args.spacing, T, U, V, P, frame=hdr, avg=args.avg)
    print("%s: %d x %d nodes, spacing %g m, %d records (%s - %s), u10 v10%s" % (args.output, nx, ny, args.spacing, len(keep),
          epoch_to_str(T[0]), epoch_to_str(T[-1]), " pa" if P is not None else ""))


def uniform(args):
    T, U, V, P = [], [], [], []
    for line in open(args.input):
        line = line.split("$")[0].split("#")[0].split()
        if not line:
            continue
        T.append(stamp_to_epoch(line[0]))
        s, d = float(line[1]), math.radians(float(line[2]))
        U.append(s * math.cos(d))
        V.append(s * math.sin(d))
        P.append(float(line[3]) if len(line) > 3 else 101325.0)
    a = lambda x: np.array(x, np.float32).reshape(-1, 1, 1)
    write_r3m(args.output, 0.0, 0.0, 1.0, 1.0, T, a(U), a(V), a(P), avg=args.avg)
    print("%s: uniform series, %d records (%s - %s)" % (args.output, len(T), epoch_to_str(T[0]), epoch_to_str(T[-1])))


def analytic_fields(X, Y, t, c):
    """u, v, pa linear in x, y and t (t in s after the first record); c: 12 coefficients."""
    u = c[0] + c[1] * X + c[2] * Y + c[3] * t
    v = c[4] + c[5] * X + c[6] * Y + c[7] * t
    p = c[8] + c[9] * X + c[10] * Y + c[11] * t
    return u, v, p


ANALYTIC = [10.0, 1.0e-4, -2.0e-4, 1.0e-4, -3.0, 2.0e-4, 1.0e-4, -5.0e-5, 101000.0, -0.01, 0.02, -0.005]


def analytic(args):
    x0, x1, y0, y1 = args.domain
    nx = int(round((x1 - x0) / args.spacing)) + 1
    ny = int(round((y1 - y0) / args.spacing)) + 1
    X, Y = np.meshgrid(x0 + args.spacing * np.arange(nx), y0 + args.spacing * np.arange(ny))
    t0 = stamp_to_epoch(args.start)
    nrec = int(round(args.hours * 3600.0 / args.step)) + 1
    T = [t0 + n * args.step for n in range(nrec)]
    U = np.empty((nrec, ny, nx), np.float32)
    V, P = np.empty_like(U), np.empty_like(U)
    for n in range(nrec):
        U[n], V[n], P[n] = analytic_fields(X, Y, n * args.step, ANALYTIC)
    write_r3m(args.output, x0, y0, args.spacing, args.spacing, T, U, V, P)
    print("%s: analytic linear field, %d x %d nodes, %d records" % (args.output, nx, ny, nrec))


def info(args):
    m = read_r3m(args.input)
    print("%s: version %d, %d x %d nodes from (%g, %g) m, spacing %g x %g m, %d records, nvar %d, moving %d"
          % (args.input, m["version"], m["nx"], m["ny"], m["x0"], m["y0"], m["dx"], m["dy"], m["nrec"], m["nvar"], m["moving"]))
    if m["frame"] == 1:
        print("frame: transverse Mercator about %g E, %g N, rotation %g deg, X0 %g, Y0 %g" % (m["lon0"], m["lat0"], m["rot"], m["X0"], m["Y0"]))
    elif m["frame"] == 2:
        print("frame: EPSG:%d, rotation %g deg, X0 %g, Y0 %g" % (m["epsg"], m["rot"], m["X0"], m["Y0"]))
    for n in range(m["nrec"]):
        s = np.hypot(m["u"][n], m["v"][n])
        line = "  %s  U10 max %6.2f mean %6.2f m/s" % (epoch_to_str(m["time"][n]), s.max(), s.mean())
        if m["pa"] is not None:
            line += "  pa min %9.1f Pa" % m["pa"][n].min()
        print(line)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    g = sub.add_parser("grid")
    g.add_argument("input")
    g.add_argument("-o", "--output", default="met-forcing.r3m")
    g.add_argument("--origin", type=float, nargs=2, metavar=("LON0", "LAT0"))
    g.add_argument("--crs")
    g.add_argument("--offset", type=float, nargs=2, metavar=("X0", "Y0"))
    g.add_argument("--rotation", type=float, default=0.0, help="model x axis, counter-clockwise from east [deg]")
    g.add_argument("--u", default="u10,uas")
    g.add_argument("--v", default="v10,vas")
    g.add_argument("--p", default="msl,psl", help="sea level pressure variable(s); '' for none")
    g.add_argument("--p-scale", type=float, default=1.0, help="factor to Pa (e.g. 100 for hPa)")
    g.add_argument("--time", default="time,valid_time")
    g.add_argument("--lon", default="longitude")
    g.add_argument("--lat", default="latitude")
    g.add_argument("--domain", type=float, nargs=4, required=True, metavar=("X0", "X1", "Y0", "Y1"))
    g.add_argument("--spacing", type=float, required=True)
    g.add_argument("--start")
    g.add_argument("--end")
    g.add_argument("--avg", type=float, default=3600.0, help="averaging period of u10 [s] (ERA5: hourly)")
    g.set_defaults(func=grid)

    u = sub.add_parser("uniform")
    u.add_argument("input")
    u.add_argument("-o", "--output", default="met-forcing.r3m")
    u.add_argument("--avg", type=float, default=0.0)
    u.set_defaults(func=uniform)

    a = sub.add_parser("analytic")
    a.add_argument("-o", "--output", default="met-forcing.r3m")
    a.add_argument("--domain", type=float, nargs=4, required=True, metavar=("X0", "X1", "Y0", "Y1"))
    a.add_argument("--spacing", type=float, required=True)
    a.add_argument("--start", default="20180101.000000")
    a.add_argument("--hours", type=float, default=6.0)
    a.add_argument("--step", type=float, default=3600.0)
    a.set_defaults(func=analytic)

    i = sub.add_parser("info")
    i.add_argument("input")
    i.set_defaults(func=info)

    args = ap.parse_args(argv)
    if args.cmd == "grid":
        # several candidate names for the time, lon and lat variables
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import seastate_forcing as sf
        nc = sf.NC(args.input)
        args.time = first_var(nc, args.time) or args.time.split(",")[0]
        if not nc.has(args.lon) and nc.has("lon"):
            args.lon, args.lat = "lon", "lat"
    args.func(args)


if __name__ == "__main__":
    main()
