#!/usr/bin/env python3
# Architect: Hans Bihs
"""
REEF3D::SEASTATE forcing converter: NetCDF spectra and wind -> REEF3D input files.

  spectra   point spectra E(time, station, frequency, direction) from WW3, NORA3 / NORAC (MET Norway)
            or similar NetCDF files -> seastate-boundary.spc, a nonstationary SWAN 2D spectral file
            in model coordinates (A 711 3)
  wind      10 m wind u10, v10 on a lon/lat (or projected) grid, e.g. ERA5 or NORA3 atmosphere
            -> seastate-wind.dat on a regular grid aligned with the model axes (A 730 2)
  point     lon/lat -> model coordinates (to check the set-up)

Model coordinates. The model's x axis points east rotated counter-clockwise by --rotation degrees; the
origin is either
  --origin LON0 LAT0            transverse Mercator on a sphere about LON0, LAT0 (no extra packages), or
  --crs EPSG:xxxx --offset X0 Y0 a projected CRS (e.g. UTM 33N = EPSG:32633, as the DIVEMesh bathymetry)
                                 minus an offset (needs pyproj).
Directions and wind vectors are turned from true north to the model frame, including the meridian
convergence of the projection at every point.

Spectra: the variable (default efth) in m^2/Hz/rad (WW3 "m2 s rad-1"; --spec-units deg for per degree);
direction convention from the direction variable's standard_name (sea_surface_wave_from_direction or
..._to_direction) or --dirconv. Any dimensions besides time, frequency and direction are stations
(several spatial dimensions are flattened). Station coordinates from --lon/--lat variables.

Times: the time variable with CF units "<seconds|minutes|hours|days> since YYYY-MM-DD[ HH:MM:SS]".
Written as YYYYMMDD.HHMMSS; set A 780 to the model start time (or 0: the first time in the file).

NetCDF backends: netCDF4 (NetCDF-3 and -4), otherwise scipy.io (NetCDF-3 only).

Examples
  python3 seastate_forcing.py spectra ww3_points.nc --origin 5.9 62.3 --rotation 20 -o seastate-boundary.spc
  python3 seastate_forcing.py wind era5_wind.nc --origin 5.9 62.3 --rotation 20 \\
          --domain 0 60000 0 40000 --spacing 2000 -o seastate-wind.dat
"""
import argparse
import datetime as dt
import math
import re
import sys

import numpy as np

R_EARTH = 6371000.0


# ------------------------------------------------------------------ NetCDF access
class NC:
    def __init__(self, path):
        try:
            import netCDF4
            self.ds = netCDF4.Dataset(path)
            self.kind = "netCDF4"
            self.ds.set_auto_mask(False)
        except ImportError:
            from scipy.io import netcdf_file
            try:
                self.ds = netcdf_file(path, "r", mmap=False)
            except Exception as e:
                sys.exit("cannot read %s with scipy.io (NetCDF-3 only); install netCDF4 for NetCDF-4 files: %s" % (path, e))
            self.kind = "scipy"

    def has(self, name):
        return name in self.ds.variables

    def dims(self, name):
        return tuple(self.ds.variables[name].dimensions)

    def attr(self, name, a, default=None):
        v = self.ds.variables[name]
        try:
            val = getattr(v, a)
        except AttributeError:
            return default
        if isinstance(val, bytes):
            val = val.decode()
        return val

    def get(self, name):
        v = self.ds.variables[name]
        a = np.array(v[:], dtype=float)
        sf = self.attr(name, "scale_factor", 1.0)
        ao = self.attr(name, "add_offset", 0.0)
        if self.kind == "scipy" and (sf != 1.0 or ao != 0.0):
            a = a * float(sf) + float(ao)
        fv = self.attr(name, "_FillValue")
        if fv is not None:
            a[a == float(fv) * (float(sf) if self.kind == "scipy" else 1.0) + (float(ao) if self.kind == "scipy" else 0.0)] = np.nan
        return a


def read_times(nc, name):
    units = nc.attr(name, "units", "")
    m = re.match(r"\s*(second|minute|hour|day)s?\s+since\s+(\d{4})-(\d{1,2})-(\d{1,2})[T ]?(\d{1,2})?:?(\d{1,2})?:?(\d{1,2}(\.\d*)?)?", units)
    if not m:
        sys.exit("time variable %s: units '%s' not understood (expected '<unit> since YYYY-MM-DD HH:MM:SS')" % (name, units))
    scale = {"second": 1.0, "minute": 60.0, "hour": 3600.0, "day": 86400.0}[m.group(1)]
    ref = dt.datetime(int(m.group(2)), int(m.group(3)), int(m.group(4)), int(m.group(5) or 0), int(m.group(6) or 0),
                      int(float(m.group(7) or 0)))
    return [ref + dt.timedelta(seconds=float(v) * scale) for v in np.atleast_1d(nc.get(name))]


def stamp(t):
    t = t + dt.timedelta(microseconds=500000)
    return "%04d%02d%02d.%02d%02d%02d" % (t.year, t.month, t.day, t.hour, t.minute, t.second)


# ------------------------------------------------------------------ projection to model coordinates
class Frame:
    def __init__(self, args):
        self.rot = math.radians(args.rotation)
        if args.crs:
            try:
                import pyproj
            except ImportError:
                sys.exit("--crs needs pyproj; use --origin LON0 LAT0 instead")
            self.tr = pyproj.Transformer.from_crs("EPSG:4326", args.crs, always_xy=True)
            self.inv = pyproj.Transformer.from_crs(args.crs, "EPSG:4326", always_xy=True)
            self.off = args.offset or [0.0, 0.0]
            self.mode = "crs"
        elif args.origin:
            self.lon0, self.lat0 = math.radians(args.origin[0]), math.radians(args.origin[1])
            self.mode = "tm"
        else:
            sys.exit("give --origin LON0 LAT0 or --crs EPSG:xxxx --offset X0 Y0")

    # projected (east, north) of the frame, before rotation
    def _fwd(self, lon, lat):
        lon = np.asarray(lon, float)
        lat = np.asarray(lat, float)
        if self.mode == "crs":
            e, n = self.tr.transform(lon, lat)
            return np.asarray(e) - self.off[0], np.asarray(n) - self.off[1]
        lam = np.radians(lon) - self.lon0
        phi = np.radians(lat)
        b = np.cos(phi) * np.sin(lam)
        e = 0.5 * R_EARTH * np.log((1.0 + b) / (1.0 - b))
        n = R_EARTH * (np.arctan2(np.tan(phi), np.cos(lam)) - self.lat0)
        return e, n

    def _inv(self, e, n):
        e = np.asarray(e, float)
        n = np.asarray(n, float)
        if self.mode == "crs":
            lon, lat = self.inv.transform(e + self.off[0], n + self.off[1])
            return np.asarray(lon), np.asarray(lat)
        d = n / R_EARTH + self.lat0
        phi = np.arcsin(np.sin(d) / np.cosh(e / R_EARTH))
        lam = self.lon0 + np.arctan2(np.sinh(e / R_EARTH), np.cos(d))
        return np.degrees(lam), np.degrees(phi)

    def to_model(self, lon, lat):
        e, n = self._fwd(lon, lat)
        c, s = math.cos(self.rot), math.sin(self.rot)
        return c * e + s * n, -s * e + c * n

    def to_geo(self, x, y):
        c, s = math.cos(self.rot), math.sin(self.rot)
        e = c * np.asarray(x) - s * np.asarray(y)
        n = s * np.asarray(x) + c * np.asarray(y)
        return self._inv(e, n)

    def north(self, lon, lat):
        """Direction of true north in the model frame [rad, counter-clockwise from the model x axis]."""
        lon = np.asarray(lon, float)
        lat = np.asarray(lat, float)
        x0, y0 = self.to_model(lon, lat)
        x1, y1 = self.to_model(lon, lat + 1.0e-4)
        return np.arctan2(y1 - y0, x1 - x0)


# ------------------------------------------------------------------ spectra
def spectra(args):
    nc = NC(args.input)
    fr = Frame(args)
    var = args.var
    if not nc.has(var):
        sys.exit("no variable %s in %s (use --var)" % (var, args.input))
    dims = nc.dims(var)
    E = nc.get(var)
    tdim = [d for d in dims if d == args.time]
    fdim, ddim = args.freq_dim or args.freq, args.dir_dim or args.dir
    if fdim not in dims or ddim not in dims:
        sys.exit("%s%s: frequency and direction dimensions %s, %s not found (use --freq-dim / --dir-dim)" % (var, dims, fdim, ddim))
    sdims = [d for d in dims if d not in (args.time, fdim, ddim)]
    order = ([dims.index(args.time)] if tdim else []) + [dims.index(d) for d in sdims] + [dims.index(fdim), dims.index(ddim)]
    E = np.transpose(E, order)
    if not tdim:
        E = E[None]
    nt = E.shape[0]
    E = E.reshape(nt, -1, E.shape[-2], E.shape[-1])
    ns = E.shape[1]

    f = nc.get(args.freq)
    d = nc.get(args.dir)
    lon = nc.get(args.lon)
    lat = nc.get(args.lat)
    if lon.size == ns and lat.size == ns:
        lon, lat = lon.reshape(-1), lat.reshape(-1)          # per station, or 2D lon/lat of a gridded field
    elif lon.ndim == 1 and lat.ndim == 1 and lon.size * lat.size == ns:
        LON, LAT = np.meshgrid(lon, lat)                     # 1D axes of a regular lon/lat field (lat, lon)
        lon, lat = LON.reshape(-1), LAT.reshape(-1)
    else:
        sys.exit("%d stations, but %d longitudes and %d latitudes" % (ns, lon.size, lat.size))

    st = args.stations if args.stations else list(range(ns))
    times = read_times(nc, args.time) if tdim else [dt.datetime(1970, 1, 1)]
    keep = [k for k, t in enumerate(times) if (args.start is None or stamp(t) >= args.start) and (args.end is None or stamp(t) <= args.end)]

    conv = args.dirconv
    if conv is None:
        sn = (nc.attr(args.dir, "standard_name", "") or "") + " " + (nc.attr(args.dir, "long_name", "") or "")
        if "from" in sn:
            conv = "from"
        elif "to" in sn.split() or "to_direction" in sn:
            conv = "to"
        else:
            sys.exit("direction convention of %s unknown: give --dirconv from|to" % args.dir)

    # direction of propagation in the geographic frame: counter-clockwise from east
    geo = np.radians(270.0 - d) if conv == "from" else np.radians(90.0 - d)
    unit = 1.0 if args.spec_units == "deg" else math.pi / 180.0       # -> per degree

    x, y = fr.to_model(lon[st], lat[st])
    north = fr.north(lon[st], lat[st])

    # model-frame direction per station: geo angle relative to true north, plus north in the model
    # frame; all stations get the directions of the first station (convergence differences are
    # applied by turning each spectrum to that common set by linear periodic interpolation)
    th0 = np.degrees(geo + (north[0] - math.pi / 2.0)) % 360.0
    o = np.argsort(th0)
    dirs = th0[o]

    out = open(args.output, "w")
    out.write("SWAN   1                                Swan standard spectral file, version\n")
    out.write("$   Data produced by tools/seastate_forcing.py from %s\n" % args.input)
    out.write("$   model frame: %s, rotation %.6g deg\n" % (("origin %.6f %.6f (transverse Mercator)" % tuple(args.origin)) if args.origin else
                                                         ("%s minus %s" % (args.crs, args.offset)), args.rotation))
    if tdim:
        out.write("TIME                                    time-dependent data\n     1                                  time coding option\n")
    out.write("LOCATIONS                               locations in x-y-space\n%6d                                  number of locations\n" % len(st))
    for a, b in zip(x, y):
        out.write("%16.4f %16.4f\n" % (a, b))
    out.write("AFREQ                                   absolute frequencies in Hz\n%6d                                  number of frequencies\n" % len(f))
    for v in f:
        out.write("%12.6f\n" % v)
    out.write("CDIR                                    spectral Cartesian directions in degr\n%6d                                  number of directions\n" % len(dirs))
    for v in dirs:
        out.write("%12.4f\n" % v)
    out.write("QUANT\n     1                                  number of quantities in table\nVaDens                                  variance densities in m2/Hz/degr\nm2/Hz/degr                              unit\n   -0.9900E+02                          exception value\n")

    for k in keep:
        if tdim:
            out.write("%s                         date and time\n" % stamp(times[k]))
        for n, s in enumerate(st):
            spec = np.nan_to_num(E[k, s], nan=0.0) * unit
            if n > 0 and abs(north[n] - north[0]) > 1.0e-9:
                # turn this station's spectrum by its convergence difference
                sh = math.degrees(north[n] - north[0])
                th = (th0 + sh) % 360.0
                spec = np.array([np.interp(dirs, np.sort(th), row[np.argsort(th)], period=360.0) for row in spec])
            else:
                spec = spec[:, o]
            spec = np.maximum(spec, 0.0)
            mx = spec.max()
            if not mx > 0.0:
                out.write("ZERO\n")
                continue
            fac = mx / 99999.0
            out.write("FACTOR\n%18.8E\n" % fac)
            for row in spec:
                out.write(" ".join("%6d" % int(round(v / fac)) for v in row) + "\n")
    out.close()
    print("%s: %d stations, %d times (%s - %s), %d frequencies, %d directions" %
          (args.output, len(st), len(keep), stamp(times[keep[0]]), stamp(times[keep[-1]]), len(f), len(dirs)))
    for n, s in enumerate(st):
        print("  station %d: lon %.5f lat %.5f -> x %.1f m, y %.1f m" % (s, lon[s], lat[s], x[n], y[n]))


# ------------------------------------------------------------------ wind
def wind(args):
    nc = NC(args.input)
    fr = Frame(args)
    u = nc.get(args.u)
    v = nc.get(args.v)
    dims = nc.dims(args.u)
    lon = nc.get(args.lon)
    lat = nc.get(args.lat)
    times = read_times(nc, args.time)
    if dims[0] != args.time:
        sys.exit("%s: the first dimension must be time (%s)" % (args.u, dims))
    keep = [k for k, t in enumerate(times) if (args.start is None or stamp(t) >= args.start) and (args.end is None or stamp(t) <= args.end)]

    x0, x1, y0, y1 = args.domain
    nx = int(round((x1 - x0) / args.spacing)) + 1
    ny = int(round((y1 - y0) / args.spacing)) + 1
    X, Y = np.meshgrid(x0 + args.spacing * np.arange(nx), y0 + args.spacing * np.arange(ny))
    LON, LAT = fr.to_geo(X, Y)
    north = fr.north(LON, LAT)            # true north in the model frame

    if lon.ndim == 1 and lat.ndim == 1:
        from scipy.interpolate import RegularGridInterpolator
        la = lat.copy()
        flip = la[0] > la[-1]
        lo = lon.copy()
        lonq = np.where(LON < lo.min(), LON + 360.0, LON) if lo.max() > 180.0 else LON

        def interp(field):
            fld = field[::-1, :] if flip else field
            g = RegularGridInterpolator((la[::-1] if flip else la, lo), fld, bounds_error=False, fill_value=None)
            return g(np.column_stack([LAT.reshape(-1), lonq.reshape(-1)])).reshape(LAT.shape)
    else:
        from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
        sx, sy = fr.to_model(lon.reshape(-1), lat.reshape(-1))
        pts = np.column_stack([sx, sy])
        from scipy.spatial import Delaunay
        tri = Delaunay(pts)

        def interp(field):
            vals = field.reshape(-1)
            a = LinearNDInterpolator(tri, vals)(X, Y)
            bad = np.isnan(a)
            if bad.any():
                a[bad] = NearestNDInterpolator(pts, vals)(X[bad], Y[bad])
            return a

    out = open(args.output, "w")
    out.write("REEF3D-SEASTATE wind field from %s (tools/seastate_forcing.py)\n" % args.input)
    out.write("%d %d\n%.6f %.6f %.6f %.6f\n" % (nx, ny, x0, y0, args.spacing, args.spacing))
    for k in keep:
        ue = interp(u[k])
        vn = interp(v[k])
        # (east, north) -> model frame: north points at angle 'north', east at north - 90 deg
        um = ue * np.sin(north) + vn * np.cos(north)
        vm = -ue * np.cos(north) + vn * np.sin(north)
        out.write("%s\n" % stamp(times[k]))
        for comp in (um, vm):
            for j in range(ny):
                out.write(" ".join("%.3f" % val for val in comp[j]) + "\n")
    out.close()
    print("%s: %d x %d nodes, spacing %g m, %d times (%s - %s)" % (args.output, nx, ny, args.spacing, len(keep), stamp(times[keep[0]]), stamp(times[keep[-1]])))


def point(args):
    fr = Frame(args)
    x, y = fr.to_model(args.lon, args.lat)
    print("lon %.6f lat %.6f -> x %.2f m, y %.2f m; true north at %.3f deg in the model frame"
          % (args.lon, args.lat, x, y, math.degrees(fr.north(args.lon, args.lat))))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def frame(p):
        p.add_argument("--origin", type=float, nargs=2, metavar=("LON0", "LAT0"))
        p.add_argument("--crs")
        p.add_argument("--offset", type=float, nargs=2, metavar=("X0", "Y0"))
        p.add_argument("--rotation", type=float, default=0.0, help="model x axis, counter-clockwise from east [deg]")

    s = sub.add_parser("spectra")
    s.add_argument("input")
    s.add_argument("-o", "--output", default="seastate-boundary.spc")
    frame(s)
    s.add_argument("--var", default="efth")
    s.add_argument("--freq", default="frequency")
    s.add_argument("--dir", default="direction")
    s.add_argument("--freq-dim")
    s.add_argument("--dir-dim")
    s.add_argument("--time", default="time")
    s.add_argument("--lon", default="longitude")
    s.add_argument("--lat", default="latitude")
    s.add_argument("--stations", type=int, nargs="*")
    s.add_argument("--dirconv", choices=["from", "to"])
    s.add_argument("--spec-units", choices=["rad", "deg"], default="rad")
    s.add_argument("--start")
    s.add_argument("--end")
    s.set_defaults(func=spectra)

    w = sub.add_parser("wind")
    w.add_argument("input")
    w.add_argument("-o", "--output", default="seastate-wind.dat")
    frame(w)
    w.add_argument("--u", default="u10")
    w.add_argument("--v", default="v10")
    w.add_argument("--time", default="time")
    w.add_argument("--lon", default="longitude")
    w.add_argument("--lat", default="latitude")
    w.add_argument("--domain", type=float, nargs=4, required=True, metavar=("X0", "X1", "Y0", "Y1"))
    w.add_argument("--spacing", type=float, required=True)
    w.add_argument("--start")
    w.add_argument("--end")
    w.set_defaults(func=wind)

    q = sub.add_parser("point")
    q.add_argument("lon", type=float)
    q.add_argument("lat", type=float)
    frame(q)
    q.set_defaults(func=point)

    args = ap.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
