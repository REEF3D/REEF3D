#!/usr/bin/env python3
# Architect: Hans Bihs
"""
Test data of the blend W 200 4 for met_test.cpp: a background met file (met_blend_bg.r3m), a moving nest
(met_blend_nest.r3m) about the track met_vortex_nh.trk, and the composite of tools/met_forcing.py (class Composite,
the reference of src/met_blend.cpp) at 120 points and 4 times for four parameter sets (met_blend_ref.txt).

  python3 met_blend_make.py        (in Regression/unit)
"""
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools"))
import met_forcing as mf
import met_track as mt

FRAME = (-86.0, 26.0, 10.0, 5000.0, -3000.0)          # W 205 lon0 lat0 rot, W 206 X0 Y0
TRACK = os.path.join(HERE, "met_vortex_nh.trk")
fr = mf.GeoFrame(*FRAME)
hdr = dict(kind=1, lon0=FRAME[0], lat0=FRAME[1], rot=FRAME[2], X0=FRAME[3], Y0=FRAME[4])
recs, rho = mt.read_trk(TRACK)


def centre_xy(t):
    lon, lat = mt.state(recs, t)[:2]
    return fr.to_model(lon, lat)


# background: 41 x 41 nodes, 40 km, 6-hourly, a smooth flow and a broad low that follows the storm
x0, y0, dx = -800e3, -600e3, 40e3
X, Y = np.meshgrid(x0 + dx * np.arange(41), y0 + dx * np.arange(41))
T = [mf.stamp_to_epoch(s) for s in ("20181009.180000", "20181010.000000", "20181010.060000", "20181010.120000",
                                     "20181010.180000")]
U, V, P = [], [], []
for n, t in enumerate(T):
    xc, yc = centre_xy(t)
    U.append(4.0 + 3.0 * np.sin(X / 300e3) + 0.5 * n)
    V.append(-2.0 + 2.0 * np.cos(Y / 250e3) - 0.3 * n)
    P.append(101000.0 + 400.0 * np.sin((X + Y) / 500e3) - 1500.0 * np.exp(-((X - xc) ** 2 + (Y - yc) ** 2) / 300e3 ** 2) + 20.0 * n)
mf.write_r3m(os.path.join(HERE, "met_blend_bg.r3m"), x0, y0, dx, dx, T, U, V, P, frame=hdr, avg=3600.0)

# nest: 21 x 21 nodes, 10 km, 3-hourly, origin on the storm (multiples of 10 km), a sharper field
TN = [mf.stamp_to_epoch("20181010.000000") + 3 * 3600.0 * k for k in range(7)]
I, J = np.meshgrid(10e3 * np.arange(21), 10e3 * np.arange(21))
U, V, P, ORG = [], [], [], []
for k, t in enumerate(TN):
    xc, yc = centre_xy(t)
    ox, oy = math.floor((xc - 100e3) / 10e3) * 10e3, math.floor((yc - 100e3) / 10e3) * 10e3
    XX, YY = ox + I, oy + J
    r2 = (XX - xc) ** 2 + (YY - yc) ** 2
    U.append(-25.0 * (YY - yc) / 60e3 * np.exp(-r2 / 80e3 ** 2) + 2.0)
    V.append(25.0 * (XX - xc) / 60e3 * np.exp(-r2 / 80e3 ** 2) - 1.0)
    P.append(100800.0 - 3000.0 * np.exp(-r2 / 60e3 ** 2) + 10.0 * k)
    ORG.append((ox, oy))
mf.write_r3m(os.path.join(HERE, "met_blend_nest.r3m"), ORG[0][0], ORG[0][1], 10e3, 10e3, TN, U, V, P, frame=hdr, origins=ORG)

# points: half near the track, half anywhere on the background
rng = np.random.default_rng(20181010)
pts = []
for k in range(60):
    t = T[1] + (T[-1] - T[1]) * rng.random()
    xc, yc = centre_xy(t)
    a, d = 2 * math.pi * rng.random(), 300e3 * rng.random()
    pts.append((round(xc + d * math.cos(a), 3), round(yc + d * math.sin(a), 3)))
for k in range(60):
    pts.append((round(-750e3 + 1500e3 * rng.random(), 3), round(-550e3 + 1500e3 * rng.random(), 3)))
times = [mf.stamp_to_epoch(s) for s in ("20181009.200000", "20181010.030000", "20181010.091500", "20181010.233000")]

VARIANTS = {
    "A": dict(nest=1, vortex=1),
    "B": dict(nest=1, vortex=1, shape=2, pressure=0, basis=2, r1=0.5, r2=1.0, vortex_time=0.0, nest_width=30e3,
              nest_time=0.0, model=1, trans=2, nest_relative=1),
    "C": dict(nest=0, vortex=1, basis=3, r1=80e3, r2=250e3, model=3, inflow=0, env=1),
    "D": dict(nest=1, vortex=0, nest_relative=1, nest_width=20e3),
}
with open(os.path.join(HERE, "met_blend_ref.txt"), "w") as out:
    out.write("# variant t x y u10 v10 pa wn wv  (Regression/unit/met_blend_make.py, tools/met_forcing.py Composite)\n")
    for name, par in VARIANTS.items():
        c = mf.Composite(fr, os.path.join(HERE, "met_blend_bg.r3m"), os.path.join(HERE, "met_blend_nest.r3m"), TRACK, **par)
        for t in times:
            c.advance(t)
            for x, y in pts:
                u, v, pa = c.at(x, y)
                wn, wv = c.weights(x, y)
                out.write("%s %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n" % (name, t, x, y, u, v, pa, wn, wv))
print("met_blend_bg.r3m, met_blend_nest.r3m, met_blend_ref.txt: %d points x %d times x %d variants" % (len(pts), len(times), len(VARIANTS)))
