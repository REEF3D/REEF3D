#!/usr/bin/env python3
# Architect: Hans Bihs
"""
geo.dat (x y z bed points) for the Berkhoff, Booy & Radder (1982) elliptic shoal.

Bathymetry as in Basilisk examples/shoal.c (refdata/waves/berkhoff_shoal/geometry.txt):
  x' = x cos20 - y sin20,  y' = x sin20 + y cos20     (shoal frame, origin at the shoal centre)
  h  = 0.45                         for x' < -5.82
  h  = 0.45 - 0.02 (5.82 + x')      otherwise (1:50 slope)
  minus the shoal  -0.3 + 0.5 sqrt(1 - (x'/3.75)^2 - (y'/5)^2)  inside (x'/3)^2 + (y'/4)^2 <= 1
  h >= 0.07 m (minimum depth as in NHWAVE/FUNWAVE set-ups, keeps the far end wet)
The REEF3D domain is the shoal frame shifted by (+10, +10): x = 0..25 m, y = 0..20 m, still water
level z = 0.45 m, so the bed is at z = 0.45 - h.

usage: berkhoff_geo.py [spacing]    (writes ./geo.dat)
"""
import math
import sys

h0, swl = 0.45, 0.45
ds = float(sys.argv[1]) if len(sys.argv) > 1 else 0.05
ca, sa = math.cos(math.radians(20.0)), math.sin(math.radians(20.0))
nx, ny = int(round(25.0 / ds)), int(round(20.0 / ds))
with open("geo.dat", "w") as f:
    for i in range(-2, nx + 3):
        for j in range(-2, ny + 3):
            X, Y = i * ds, j * ds
            x, y = X - 10.0, Y - 10.0
            xr, yr = x * ca - y * sa, x * sa + y * ca
            z0 = (5.82 + xr) / 50.0 if xr >= -5.82 else 0.0
            zs = 0.0
            if (xr / 3.0) ** 2 + (yr / 4.0) ** 2 <= 1.0:
                zs = -0.3 + 0.5 * math.sqrt(1.0 - (xr / 3.75) ** 2 - (yr / 5.0) ** 2)
            h = max(h0 - z0 - zs, 0.07)
            f.write("%.4f %.4f %.5f\n" % (X, Y, swl - h))
