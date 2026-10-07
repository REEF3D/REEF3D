#!/usr/bin/env python3
# Architect: Hans Bihs
"""
geo.dat for Thacker's parabolic bowl (SWASHES 4.2.1): z = h0 (x - L/2)^2 / a^2 with a = 1 m,
h0 = 0.5 m, L = 4 m (bed relative to the bowl bottom z = 0), 1D in x (three rows in y).

usage: thacker_geo.py [ds]
"""
import sys

ds = float(sys.argv[1]) if len(sys.argv) > 1 else 0.005
n = int(round(4.0 / ds))
with open("geo.dat", "w") as f:
    for i in range(-2, n + 3):
        x = i * ds
        for y in (-0.05, 0.0, 0.05, 0.1):
            f.write("%.5f %.3f %.6f\n" % (x, y, 0.5 * (x - 2.0) ** 2))
