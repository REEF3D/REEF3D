#!/usr/bin/env python3
# Architect: Hans Bihs
"""
REEF3D atmospheric forcing (W 200 2, 4): writer, reader and converters of the binary met file met-forcing.r3m.

  grid      10 m wind u10, v10 and mean sea level pressure on a lon/lat (or projected) grid from NetCDF, e.g. ERA5
            (u10, v10, msl), NORA3 or BARRA-C2 (uas, vas, psl) -> met-forcing.r3m on a regular grid in model
            coordinates; wind vectors turned from (east, north) to the model axes incl. the meridian convergence
  uniform   spatially uniform time series (nx = ny = 1) from a table: time U10 dir [pa]
  analytic  test fields: u10, v10, pa linear in x, y and t (exact under bilinear and linear interpolation)
  low       air pressure only: a Gaussian low on a background gradient, static or moving (inverse barometer tests)
  nest      a storm-following nest (e.g. HAFS) from NetCDF files on lon/lat grids that move from file to file ->
            met-nest.r3m, a moving grid of fixed size in model coordinates (W 200 4, W 220 1 x)
  blend     the composite of W 200 4 (background, nest, vortex; the reference of src/met_blend.cpp) at points or on
            a grid
  info      header and records of a met file

Model coordinates: as tools/seastate_forcing.py, --origin LON0 LAT0 (transverse Mercator on a sphere, the frame of
W 205 LON0 LAT0 ROT, with --w206 X0 Y0 the model coordinates of LON0 LAT0 as W 206) or --crs EPSG:xxxx --offset X0 Y0
(pyproj), and --rotation ROT of the model x axis counter-clockwise from east [deg]. The frame is stored in the file header and checked against W 205 / W 206.

File format version 1 (little-endian, src/met_file.h)
  header 144 bytes: char[8] "REEF3DMT"; int32 version, nx, ny, nvar, moving, nrec, frame, epsg;
                    float64 x0, y0, dx, dy, lon0, lat0, rot, X0, Y0, avg, reserved[3]
  records:          float64 time [s since 1970 UTC], xr0, yr0, reserved; float32 u10[ny][nx], v10[ny][nx], (pa[ny][nx])

Examples
  python3 met_forcing.py grid era5_2018_michael.nc --origin -85.5 29.5 --domain 0 400000 -200000 200000 \\
          --spacing 5000 -o met-forcing.r3m
  python3 met_forcing.py uniform series.txt -o met-forcing.r3m      # lines: YYYYMMDD.HHMMSS U10 dir [pa]
  python3 met_forcing.py analytic --domain 0 10000 0 5000 --spacing 500 --start 20180101.000000 --hours 6 --step 3600
  python3 met_forcing.py nest 'hafs.*.nc' --origin -85.5 29.5 --size 600000 600000 --spacing 3000 -o met-nest.r3m
  python3 met_forcing.py blend met-forcing.r3m --track met-vortex.trk --domain 0 400000 -200000 200000 --spacing 5000
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
        w = getattr(args, "w206", None) or (0.0, 0.0)
        hdr = dict(kind=1, lon0=args.origin[0], lat0=args.origin[1], rot=args.rotation, X0=w[0], Y0=w[1])
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
    w206 = args.w206 or (0.0, 0.0)
    LON, LAT = fr.to_geo(X - w206[0], Y - w206[1])
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
        pts = np.column_stack([sx + w206[0], sy + w206[1]])
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


def low(args):
    """Air pressure only (no wind): a Gaussian low pa = pref + gx x + gy y - dp exp(-r^2/R^2), its centre
    moving from (xc, yc) with (vx, vy); for the air pressure gradient and the inverse barometer (W 208, W 210)."""
    x0, x1, y0, y1 = args.domain
    nx = int(round((x1 - x0) / args.spacing)) + 1
    ny = int(round((y1 - y0) / args.spacing)) + 1
    X, Y = np.meshgrid(x0 + args.spacing * np.arange(nx), y0 + args.spacing * np.arange(ny))
    t0 = stamp_to_epoch(args.start)
    nrec = int(round(args.hours * 3600.0 / args.step)) + 1
    T = [t0 + n * args.step for n in range(nrec)]
    U = np.zeros((nrec, ny, nx), np.float32)
    P = np.empty_like(U)
    for n in range(nrec):
        t = n * args.step
        xc, yc = args.centre[0] + args.velocity[0] * t, args.centre[1] + args.velocity[1] * t
        r2 = (X - xc) ** 2 + (Y - yc) ** 2
        P[n] = args.pref + args.gradient[0] * X + args.gradient[1] * Y - args.dp * np.exp(-r2 / args.radius ** 2)
    write_r3m(args.output, x0, y0, args.spacing, args.spacing, T, U, U, P)
    print("%s: pressure low, %d x %d nodes, %d records, pa %.1f - %.1f Pa" % (args.output, nx, ny, nrec, P.min(), P.max()))


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



# ------------------------------------------------------------------ reference of the C++ readers and the blend (W 200 4)
class GeoFrame:
    """met_frame of REEF3D operation by operation (W 205 lon0 lat0 rot, W 206 X0 Y0), scalar"""
    R = 6371000.0

    def __init__(self, lon0, lat0, rot=0.0, X0=0.0, Y0=0.0):
        rad = math.pi / 180.0
        self.lon0, self.lat0, self.rot, self.X0, self.Y0 = lon0, lat0, rot, X0, Y0
        self.lam0, self.phi0 = lon0 * rad, lat0 * rad
        self.c, self.s = math.cos(rot * rad), math.sin(rot * rad)

    def to_model(self, lon, lat):
        rad = math.pi / 180.0
        lam = lon * rad - self.lam0
        phi = lat * rad
        b = math.cos(phi) * math.sin(lam)
        e = 0.5 * self.R * math.log((1.0 + b) / (1.0 - b))
        n = self.R * (math.atan2(math.tan(phi), math.cos(lam)) - self.phi0)
        return self.c * e + self.s * n + self.X0, -self.s * e + self.c * n + self.Y0

    def to_geo(self, x, y):
        rad = math.pi / 180.0
        xx, yy = x - self.X0, y - self.Y0
        e = self.c * xx - self.s * yy
        n = self.s * xx + self.c * yy
        d = n / self.R + self.phi0
        return (self.lam0 + math.atan2(math.sinh(e / self.R), math.cos(d))) / rad, math.asin(math.sin(d) / math.cosh(e / self.R)) / rad

    def north_geo(self, lon, lat):
        x0, y0 = self.to_model(lon, lat)
        x1, y1 = self.to_model(lon, lat + 1.0e-4)
        return math.atan2(y1 - y0, x1 - x0)


class R3M:
    """met_file of REEF3D on the whole grid: bilinear as nested lerps, linear in time, clamped"""

    def __init__(self, path):
        self.m = read_r3m(path)
        m = self.m
        self.nx, self.ny, self.nrec, self.nvar = m["nx"], m["ny"], m["nrec"], m["nvar"]
        self.dx, self.dy = m["dx"], m["dy"]
        self.times = [float(t) for t in m["time"]]
        self.k = 0

    def advance(self, t):
        k = 0
        if self.nrec > 1:
            k = int(np.searchsorted(self.times, t, side="right")) - 1
            k = min(max(k, 0), self.nrec - 2)
        self.k = k

    def weight(self, t):
        if self.nrec < 2:
            return 0.0
        t0, t1 = self.times[self.k], self.times[self.k + 1]
        return min(max((t - t0) / (t1 - t0), 0.0), 1.0)

    def _space(self, n, x, y):
        m = self.m
        xr, yr = float(m["origin"][n][0]), float(m["origin"][n][1])
        i = j = 0
        wx = wy = 0.0
        if self.nx > 1:
            fx = min(max((x - xr) / self.dx, 0.0), float(self.nx - 1))
            i = min(int(math.floor(fx)), self.nx - 2)
            wx = fx - float(i)
        if self.ny > 1:
            fy = min(max((y - yr) / self.dy, 0.0), float(self.ny - 1))
            j = min(int(math.floor(fy)), self.ny - 2)
            wy = fy - float(j)
        i1 = i + 1 if wx > 0.0 else i
        j1 = j + 1 if wy > 0.0 else j
        lerp = lambda a, b, w: a + w * (b - a)
        out = []
        for f in (m["u"], m["v"], m["pa"]):
            if f is None:
                out.append(101325.0)
                continue
            g = f[n]
            out.append(lerp(lerp(float(g[j, i]), float(g[j, i1]), wx), lerp(float(g[j1, i]), float(g[j1, i1]), wx), wy))
        return out

    def at(self, t, x, y):
        u, v, p = self._space(self.k, x, y)
        w = self.weight(t)
        if w > 0.0:
            u1, v1, p1 = self._space(self.k + 1, x, y)
            u, v, p = u + w * (u1 - u), v + w * (v1 - v), p + w * (p1 - p)
        return u, v, p

    def origin(self, t):
        o = self.m["origin"]
        ox, oy = float(o[self.k][0]), float(o[self.k][1])
        w = self.weight(t)
        if w > 0.0:
            ox, oy = ox + w * (float(o[self.k + 1][0]) - ox), oy + w * (float(o[self.k + 1][1]) - oy)
        return ox, oy

    def rec_index(self, q):
        return min(self.k + q, self.nrec - 1)

    def at_record(self, q, x, y):
        """values of the loaded record q (0, 1) at (x, y) (met_file::at_record)"""
        return self._space(self.rec_index(q), x, y)

    def sample(self, n, var, x, y):
        """one variable of record n, bilinear on the whole grid (met_file::sample)"""
        return self._space(n, x, y)[var]


def taper(s, shape):
    if s <= 0.0:
        return 0.0
    if s >= 1.0:
        return 1.0
    return s if shape == 2 else 0.5 * (1.0 - math.cos(math.pi * s))


def taper_time(t, t0, t1, T, shape):
    if t < t0 or t > t1:
        return 0.0
    if not T > 0.0:
        return 1.0
    return taper(min(t - t0, t1 - t) / T, shape)


def mix(f0, f1, w):
    if w <= 0.0:
        return f0
    if w >= 1.0:
        return f1
    return f0 + w * (f1 - f0)


class Composite:
    """reference of met_source_blend (W 200 4): background, moving nest and vortex with tapered weights.
    par: dict(nest, vortex, nest_width, nest_time, basis, r1, r2, pressure, vortex_time, shape, env,
              model, wrf, one2ten, inflow, trans) as W 220-226, W 215-219"""

    DEFAULT = dict(nest=0, vortex=1, nest_width=50000.0, nest_time=3600.0, basis=1, r1=3.0, r2=6.0, pressure=1,
                   vortex_time=21600.0, shape=1, env=0, nest_relative=0, model=2, wrf=0.9, one2ten=0.8928, inflow=1, trans=1)

    def __init__(self, frame, background, nestfile=None, trackfile=None, **par):
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import met_track
        self.mt = met_track
        self.p = dict(self.DEFAULT)
        self.p.update(par)
        self.frame = frame
        self.bg = R3M(background)
        self.nest = R3M(nestfile) if (self.p["nest"] == 1) else None
        need = self.p["vortex"] == 1 or (self.p["nest"] == 1 and self.p["nest_relative"] == 1)
        self.recs, self.rho = met_track.read_trk(trackfile) if need else (None, None)
        self.ring_cache = {}

    def radii(self, t):
        lon, lat, pc, pn, rmw, r34 = self.mt.state(self.recs, t)
        sc = rmw if self.p["basis"] == 1 else (r34 if self.p["basis"] == 2 else 1.0)
        return self.p["r1"] * sc, self.p["r2"] * sc

    def ring(self, n):
        """the background and the vortex on the ring r2 at background record n (met_source_blend::ring)"""
        if n in self.ring_cache:
            return self.ring_cache[n]
        p = self.p
        t = self.bg.times[n]
        lon, lat, pc, pn, rmw, r34 = self.mt.state(self.recs, t)
        a, b = self.radii(t)
        xc, yc = self.frame.to_model(lon, lat)
        X, Y, R = [], [], []
        su = sv = sp = spv = sx = sy = sr = 0.0
        for m in range(16):
            lon2, lat2 = self.mt.destination(lon, lat, 360.0 * float(m) / 16.0, b)
            x, y = self.frame.to_model(lon2, lat2)
            pb = self.bg.sample(n, 2, x, y)
            ub = vb = 0.0
            if p["env"] == 1:
                ub, vb = self.bg.sample(n, 0, x, y), self.bg.sample(n, 1, x, y)
            ue, vn, pv = self.mt.at_time(self.recs, t, lon2, lat2, p["model"], self.rho, wrf=p["wrf"], one2ten=p["one2ten"],
                                         inflow=p["inflow"], trans=p["trans"])
            chi = self.frame.north_geo(lon2, lat2)
            uv = ue * math.sin(chi) + vn * math.cos(chi)
            vv = -ue * math.cos(chi) + vn * math.sin(chi)
            X.append(x - xc); Y.append(y - yc); R.append(pb - pv)
            sp += pb; spv += pv; sr += R[-1]; sx += X[-1]; sy += Y[-1]
            su += ub - uv; sv += vb - vv
        rg = dict(p=sp / 16.0, pv=spv / 16.0, u=su / 16.0 if p["env"] == 1 else 0.0, v=sv / 16.0 if p["env"] == 1 else 0.0)
        mx, my, mr = sx / 16.0, sy / 16.0, sr / 16.0
        sxx = syy = sxy = sxp = syp = 0.0
        for m in range(16):
            dx, dy, dr = X[m] - mx, Y[m] - my, R[m] - mr
            sxx += dx * dx; syy += dy * dy; sxy += dx * dy; sxp += dx * dr; syp += dy * dr
        det = sxx * syy - sxy * sxy
        rg["gx"] = (sxp * syy - syp * sxy) / det if det > 0.0 else 0.0
        rg["gy"] = (syp * sxx - sxp * sxy) / det if det > 0.0 else 0.0
        self.ring_cache[n] = rg
        return rg

    def advance(self, t):
        p = self.p
        self.t = t
        self.bg.advance(t)
        self.wn_t = 0.0
        if self.nest is not None:
            self.nest.advance(t)
            self.wn_t = taper_time(t, self.nest.times[0], self.nest.times[-1], p["nest_time"], p["shape"])
            self.nx0, self.ny0 = self.nest.origin(t)
            self.nlx, self.nly = float(self.nest.nx - 1) * self.nest.dx, float(self.nest.ny - 1) * self.nest.dy
            if p["nest_relative"] == 1:
                xc, yc = self.frame.to_model(*self.mt.state(self.recs, t)[:2])
                self.nest_w = self.nest.weight(t)
                self.nest_d, o = [], []
                for q in range(2):
                    n = self.nest.rec_index(q)
                    xq, yq = self.frame.to_model(*self.mt.state(self.recs, self.nest.times[n])[:2])
                    d = (xc - xq, yc - yq)
                    self.nest_d.append(d)
                    org = self.nest.m["origin"][n]
                    o.append((float(org[0]) + d[0], float(org[1]) + d[1]))
                w = self.nest_w
                self.nx0 = o[0][0] + w * (o[1][0] - o[0][0]) if w > 0.0 else o[0][0]
                self.ny0 = o[0][1] + w * (o[1][1] - o[0][1]) if w > 0.0 else o[0][1]
        self.wv_t = 0.0
        if self.recs is not None and p["vortex"] == 1:
            self.wv_t = taper_time(t, self.recs[0]["t"], self.recs[-1]["t"], p["vortex_time"], p["shape"])
            self.lonc, self.latc, self.pc, self.pn, self.rmw, self.r34 = self.mt.state(self.recs, t)
            self.r1m, self.r2m = self.radii(t)
            self.xcm, self.ycm = self.frame.to_model(self.lonc, self.latc)
            self.env = dict(p=self.pn, pv=self.pn, u=0.0, v=0.0, gx=0.0, gy=0.0)
            if p["pressure"] == 1 or p["env"] == 1:
                k, w = 0, 0.0
                if self.bg.nrec > 1:
                    k = self.bg.k
                    w = self.bg.weight(t)
                k1 = k + 1 if self.bg.nrec > 1 else k
                a, b = self.ring(k), self.ring(k1)
                self.env = {q: (a[q] + w * (b[q] - a[q]) if w > 0.0 else a[q]) for q in ("p", "pv", "u", "v", "gx", "gy")}
            self.pb_ring = self.env["p"]

    def weights(self, x, y):
        """nest and vortex weights at (x, y) and the time of the last advance"""
        wn = wv = 0.0
        if self.wn_t > 0.0:
            d = min(min(x - self.nx0, self.nx0 + self.nlx - x), min(y - self.ny0, self.ny0 + self.nly - y))
            wn = self.wn_t * taper(d / self.p["nest_width"], self.p["shape"])
        if self.wv_t > 0.0:
            lon, lat = self.frame.to_geo(x, y)
            r, _ = self.mt.dist_bearing(self.lonc, self.latc, lon, lat)
            wv = self.wv_t * taper((self.r2m - r) / (self.r2m - self.r1m), self.p["shape"])
        return wn, wv

    def at(self, x, y):
        p, t = self.p, self.t
        u, v, pa = self.bg.at(t, x, y)
        wn, wv = self.weights(x, y)
        if wn > 0.0:
            if p["nest_relative"] == 1:
                d0, d1 = self.nest_d
                un, vn, pn = self.nest.at_record(0, x - d0[0], y - d0[1])
                if self.nest_w > 0.0:
                    u1, v1, p1 = self.nest.at_record(1, x - d1[0], y - d1[1])
                    w = self.nest_w
                    un, vn, pn = un + w * (u1 - un), vn + w * (v1 - vn), pn + w * (p1 - pn)
            else:
                un, vn, pn = self.nest.at(t, x, y)
            u, v = mix(u, un, wn), mix(v, vn, wn)
            if self.nest.nvar == 3:
                pa = mix(pa, pn, wn)
        if wv > 0.0:
            lon, lat = self.frame.to_geo(x, y)
            ue, vn, pv = self.mt.at_time(self.recs, t, lon, lat, p["model"], self.rho, wrf=p["wrf"], one2ten=p["one2ten"],
                                         inflow=p["inflow"], trans=p["trans"])
            chi = self.frame.north_geo(lon, lat)
            uv = ue * math.sin(chi) + vn * math.cos(chi)
            vv = -ue * math.cos(chi) + vn * math.sin(chi)
            if p["pressure"] == 1:
                pv = self.pc + (pv - self.pc) * max((self.env["p"] - self.pc) / (self.env["pv"] - self.pc), 0.0)
            uw, vw = uv, vv
            if p["env"] == 1:
                uw, vw = uv + self.env["u"], vv + self.env["v"]
                pv = pv + self.env["gx"] * (x - self.xcm) + self.env["gy"] * (y - self.ycm)
            u, v, pa = mix(u, uw, wv), mix(v, vw, wv), mix(pa, pv, wv)
        return u, v, pa


def blend_par(args):
    return dict(nest=1 if args.nest else 0, vortex=1 if args.track else 0, nest_width=args.nest_taper[0],
                nest_time=args.nest_taper[1], basis=int(args.radii[0]), r1=args.radii[1], r2=args.radii[2],
                pressure=args.pressure, vortex_time=args.vortex_taper, shape=args.shape, env=args.env, nest_relative=args.nest_relative, model=args.model,
                wrf=args.wrf, one2ten=args.one2ten, inflow=args.inflow, trans=args.trans)


def blend(args):
    """the composite of W 200 4 at points (reference table) or on a grid (met file, e.g. for plots or SEASTATE)"""
    m = read_r3m(args.background)
    if m["frame"] != 1 and args.frame is None:
        sys.exit("%s has no transverse Mercator frame; give --frame LON0 LAT0 ROT X0 Y0" % args.background)
    fr = GeoFrame(*args.frame) if args.frame else GeoFrame(m["lon0"], m["lat0"], m["rot"], m["X0"], m["Y0"])
    c = Composite(fr, args.background, args.nest, args.track, **blend_par(args))
    if args.points:
        pts = [tuple(map(float, l.split()[:2])) for l in open(args.points) if l.strip() and not l.startswith("#")]
        times = [stamp_to_epoch(s) for s in args.time]
        out = open(args.output, "w") if args.output else sys.stdout
        out.write("# t x y u10 v10 pa wn wv (tools/met_forcing.py blend)\n")
        for t in times:
            c.advance(t)
            for x, y in pts:
                u, v, pa = c.at(x, y)
                wn, wv = c.weights(x, y)
                out.write("%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n" % (t, x, y, u, v, pa, wn, wv))
        return
    x0, x1, y0, y1 = args.domain
    nx = int(round((x1 - x0) / args.spacing)) + 1
    ny = int(round((y1 - y0) / args.spacing)) + 1
    t0 = stamp_to_epoch(args.start) if args.start else c.bg.times[0]
    t1 = stamp_to_epoch(args.end) if args.end else c.bg.times[-1]
    T = list(np.arange(t0, t1 + 0.5 * args.step, args.step))
    U = np.empty((len(T), ny, nx), np.float32)
    V, P = np.empty_like(U), np.empty_like(U)
    for n, t in enumerate(T):
        c.advance(t)
        for j in range(ny):
            for i in range(nx):
                U[n, j, i], V[n, j, i], P[n, j, i] = c.at(x0 + i * args.spacing, y0 + j * args.spacing)
    hdr = dict(kind=1, lon0=fr.lon0, lat0=fr.lat0, rot=fr.rot, X0=fr.X0, Y0=fr.Y0)
    write_r3m(args.output or "met-blend.r3m", x0, y0, args.spacing, args.spacing, T, U, V, P, frame=hdr, avg=m["avg"])
    print("%s: composite, %d x %d nodes, %d records, U10 max %.1f m/s, pa min %.1f hPa"
          % (args.output or "met-blend.r3m", nx, ny, len(T), float(np.hypot(U, V).max()), float(P.min()) / 100.0))


def nest(args):
    """storm-following nest (e.g. HAFS) from NetCDF files -> a moving met file met-nest.r3m: per time a grid of
    fixed size in model coordinates about the centre of the source grid (or of the vortex track), origin on
    multiples of the spacing"""
    sf, fr, hdr = frame_of(args)
    files = sorted(sum([sorted(__import__("glob").glob(g)) for g in args.input], []))
    if not files:
        sys.exit("no input files")
    recs = None
    if args.track:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import met_track
        recs, _ = met_track.read_trk(args.track)
    T, U, V, P, ORG = [], [], [], [], []
    w206 = args.w206 or (0.0, 0.0)
    lx, ly = args.size
    nx = int(round(lx / args.spacing)) + 1
    ny = int(round(ly / args.spacing)) + 1
    I, J = np.meshgrid(np.arange(nx) * args.spacing, np.arange(ny) * args.spacing)
    from scipy.interpolate import RegularGridInterpolator
    for path in files:
        nc = sf.NC(path)
        tname = first_var(nc, args.time) or args.time.split(",")[0]
        lonn = args.lon if nc.has(args.lon) else ("lon" if nc.has("lon") else args.lon)
        latn = args.lat if nc.has(args.lat) else ("lat" if nc.has("lat") else args.lat)
        times = sf.read_times(nc, tname)
        vu, vv = first_var(nc, args.u), first_var(nc, args.v)
        vp = first_var(nc, args.p) if args.p else None
        if vu is None or vv is None:
            sys.exit("%s: no wind variables %s / %s" % (path, args.u, args.v))
        u, v = nc.get(vu), nc.get(vv)
        pr = nc.get(vp) if vp else None
        lon_all, lat_all = nc.get(lonn), nc.get(latn)
        if u.ndim == 2:
            u, v = u[None], v[None]
            pr = pr[None] if pr is not None else None
        for k, tk in enumerate(times):
            lon = lon_all[k] if lon_all.ndim == 2 else lon_all     # (time, x) for moving 1D axes
            lat = lat_all[k] if lat_all.ndim == 2 else lat_all
            lon = np.where(lon > 180.0, lon - 360.0, lon)
            te = datetime_to_epoch(tk)
            if recs is not None:
                lonc, latc = met_track.state(recs, te)[:2]
            else:
                lonc, latc = 0.5 * (lon.min() + lon.max()), 0.5 * (lat.min() + lat.max())
            xc, yc = fr.to_model(lonc, latc)
            xc, yc = float(xc) + w206[0], float(yc) + w206[1]
            ox = math.floor((float(xc) - 0.5 * lx) / args.spacing) * args.spacing
            oy = math.floor((float(yc) - 0.5 * ly) / args.spacing) * args.spacing
            LON, LAT = fr.to_geo(ox + I - w206[0], oy + J - w206[1])
            if LON.min() < lon.min() or LON.max() > lon.max() or LAT.min() < lat.min() or LAT.max() > lat.max():
                sys.exit("%s, %s: the nest grid of %g x %g m about (%.2f, %.2f) leaves the source grid "
                         "(%.2f..%.2f E, %.2f..%.2f N); use a smaller --size" % (path, sf.stamp(tk), lx, ly, lonc, latc,
                                                                                 lon.min(), lon.max(), lat.min(), lat.max()))
            flip = lat[0] > lat[-1]
            la = lat[::-1] if flip else lat
            north = fr.north(LON, LAT)

            def interp(field):
                fld = field[::-1, :] if flip else field
                g = RegularGridInterpolator((la, lon), fld)
                return g(np.column_stack([LAT.reshape(-1), LON.reshape(-1)])).reshape(LAT.shape)
            ue, vn = interp(u[k]), interp(v[k])
            U.append(ue * np.sin(north) + vn * np.cos(north))
            V.append(-ue * np.cos(north) + vn * np.sin(north))
            if pr is not None:
                P.append(interp(pr[k]) * args.p_scale)
            T.append(te)
            ORG.append((ox, oy))
    order = np.argsort(T)
    T = [T[i] for i in order]
    U, V = np.array(U)[order], np.array(V)[order]
    P = np.array(P)[order] if P else None
    ORG = [ORG[i] for i in order]
    write_r3m(args.output, ORG[0][0], ORG[0][1], args.spacing, args.spacing, T, U, V, P, frame=hdr, avg=args.avg, origins=ORG)
    print("%s: moving nest, %d x %d nodes, spacing %g m, %d records (%s - %s)%s" % (args.output, nx, ny, args.spacing, len(T),
          epoch_to_str(T[0]), epoch_to_str(T[-1]), ", u10 v10 pa" if P is not None else ", u10 v10"))

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
    g.add_argument("--w206", type=float, nargs=2, metavar=("X0", "Y0"), help="model coordinates of LON0 LAT0 (W 206)")
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

    l = sub.add_parser("low")
    l.add_argument("-o", "--output", default="met-forcing.r3m")
    l.add_argument("--domain", type=float, nargs=4, required=True, metavar=("X0", "X1", "Y0", "Y1"))
    l.add_argument("--spacing", type=float, required=True)
    l.add_argument("--dp", type=float, default=3000.0, help="central pressure deficit [Pa]")
    l.add_argument("--radius", type=float, default=50000.0, help="e-folding radius R [m]")
    l.add_argument("--centre", type=float, nargs=2, default=(0.0, 0.0), metavar=("XC", "YC"))
    l.add_argument("--velocity", type=float, nargs=2, default=(0.0, 0.0), metavar=("VX", "VY"))
    l.add_argument("--gradient", type=float, nargs=2, default=(0.0, 0.0), metavar=("GX", "GY"), help="background gradient [Pa/m]")
    l.add_argument("--pref", type=float, default=101325.0)
    l.add_argument("--start", default="20181010.120000")
    l.add_argument("--hours", type=float, default=6.0)
    l.add_argument("--step", type=float, default=3600.0)
    l.set_defaults(func=low)

    b = sub.add_parser("blend")
    b.add_argument("background", help="background met file (met-forcing.r3m)")
    b.add_argument("--nest", help="nest met file (met-nest.r3m), W 220 1 x")
    b.add_argument("--track", help="vortex track (met-vortex.trk), W 220 x 1")
    b.add_argument("--frame", type=float, nargs=5, metavar=("LON0", "LAT0", "ROT", "X0", "Y0"),
                   help="W 205, W 206 (default: the frame in the header of the background)")
    b.add_argument("--nest-taper", type=float, nargs=2, default=(50000.0, 3600.0), metavar=("WIDTH", "T"), help="W 221")
    b.add_argument("--radii", type=float, nargs=3, default=(1, 3.0, 6.0), metavar=("BASIS", "R1", "R2"), help="W 222")
    b.add_argument("--pressure", type=int, default=1, help="W 223")
    b.add_argument("--vortex-taper", type=float, default=21600.0, help="W 224 [s]")
    b.add_argument("--shape", type=int, default=1, help="W 225: 1 cosine, 2 linear")
    b.add_argument("--env", type=int, default=0, help="W 226: 1 the background flow and pressure gradient of the ring added to the vortex")
    b.add_argument("--nest-relative", type=int, default=0, help="W 227: 1 storm-relative time interpolation of the nest (needs --track)")
    b.add_argument("--model", type=int, default=2, help="W 215")
    b.add_argument("--wrf", type=float, default=0.9, help="W 216")
    b.add_argument("--one2ten", type=float, default=0.8928, help="W 217")
    b.add_argument("--inflow", type=int, default=1, help="W 218")
    b.add_argument("--trans", type=int, default=1, help="W 219")
    b.add_argument("--points", help="file with lines x y [m]: write the composite at the points (--time) instead of a grid")
    b.add_argument("--time", nargs="+", default=[], help="YYYYMMDD.HHMMSS (with --points)")
    b.add_argument("--domain", type=float, nargs=4, metavar=("X0", "X1", "Y0", "Y1"))
    b.add_argument("--spacing", type=float)
    b.add_argument("--start")
    b.add_argument("--end")
    b.add_argument("--step", type=float, default=3600.0)
    b.add_argument("-o", "--output")
    b.set_defaults(func=blend)

    n = sub.add_parser("nest")
    n.add_argument("input", nargs="+", help="NetCDF files (globs) with u10, v10 (msl) on a lon/lat grid that may move from file to file")
    n.add_argument("-o", "--output", default="met-nest.r3m")
    n.add_argument("--origin", type=float, nargs=2, metavar=("LON0", "LAT0"))
    n.add_argument("--crs")
    n.add_argument("--offset", type=float, nargs=2, metavar=("X0", "Y0"))
    n.add_argument("--rotation", type=float, default=0.0)
    n.add_argument("--w206", type=float, nargs=2, metavar=("X0", "Y0"), help="model coordinates of LON0 LAT0 (W 206)")
    n.add_argument("--size", type=float, nargs=2, required=True, metavar=("LX", "LY"), help="nest grid in model coordinates [m]")
    n.add_argument("--spacing", type=float, required=True)
    n.add_argument("--track", help="centre the nest grid on the vortex of this track instead of the source grid")
    n.add_argument("--u", default="u10,UGRD_10maboveground,uas")
    n.add_argument("--v", default="v10,VGRD_10maboveground,vas")
    n.add_argument("--p", default="msl,PRMSL_meansealevel,MSLET_meansealevel,psl")
    n.add_argument("--p-scale", type=float, default=1.0)
    n.add_argument("--time", default="time,valid_time")
    n.add_argument("--lon", default="longitude")
    n.add_argument("--lat", default="latitude")
    n.add_argument("--avg", type=float, default=0.0)
    n.set_defaults(func=nest)

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
