#!/usr/bin/env python3
# Architect: Hans Bihs
"""
geo.dat for the conical island of Briggs et al. (1995) (NTHMP benchmark BP6).

Island (FUNWAVE-TVD benchmarks/car_conical_island/input/conical.m): toe radius 3.6 m, side slope 1:4,
height 0.625 m (flat top), tank depth 0.32 m. In the REEF3D domain (x 0..26 m, y 0..27.6 m) the
centre is at (15, 13.8) and the bed is z = clip((3.6 - r)/4, 0, 0.625) above the flat bottom z = 0.
Dense points (spacing ds) for r < 4 m, a 0.5 m grid elsewhere.

usage: conical_geo.py [ds]
"""
import math
import sys

ds = float(sys.argv[1]) if len(sys.argv) > 1 else 0.05
xc, yc, r1, s, hc = 15.0, 13.8, 3.6, 0.25, 0.625


def zb(x, y):
    r = math.hypot(x - xc, y - yc)
    return min(max((r1 - r) * s, 0.0), hc)


pts = {}
n = int(round(26.0 / 0.5))
m = int(round(27.6 / 0.5))
for i in range(-1, n + 2):
    for j in range(-1, m + 2):
        x, y = i * 0.5, j * 0.5
        if math.hypot(x - xc, y - yc) >= 4.0:
            pts[(round(x, 4), round(y, 4))] = 0.0
k = int(round(4.2 / ds))
for i in range(-k, k + 1):
    for j in range(-k, k + 1):
        x, y = xc + i * ds, yc + j * ds
        if math.hypot(x - xc, y - yc) < 4.2:
            pts[(round(x, 4), round(y, 4))] = zb(x, y)
with open("geo.dat", "w") as f:
    for (x, y), z in pts.items():
        f.write("%.4f %.4f %.5f\n" % (x, y, z))
