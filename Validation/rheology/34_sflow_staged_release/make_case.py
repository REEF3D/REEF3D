#!/usr/bin/env python3
# 34: a block of Coulomb material on a 20 deg plane, released in two halves (F 74): the toe half (x 50-70 m) at the
# start, the back half (x 30-50 m) at t = 10 s. Writes geo.dat, sflow-level.dat and sflow-release.dat into
# cases/staged and cases/instant (instant: the same material without F 74).
import os
import numpy as np

dx, L, W = 2.0, 200.0, 40.0
x = np.arange(dx/2, L, dx); y = np.arange(dx/2, W, dx)
X, Y = np.meshgrid(x, y)
bed = np.tan(np.radians(20.0))*(L + 20.0 - X)
inblock = (X > 30) & (X < 70) & (Y > 10) & (Y < 30)
level = np.where(inblock, bed + 3.0, -9999.0)
trel = np.where(inblock & (X < 50), 10.0, -9999.0)

for c in ("staged", "instant"):
    d = os.path.join(os.path.dirname(os.path.abspath(__file__)), "cases", c)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "geo.dat"), "w") as f:
        for j in range(len(y)):
            for i in range(len(x)):
                f.write("%.2f %.2f %.3f\n" % (x[i], y[j], bed[j, i]))
    for name, v in (("sflow-level.dat", level), ("sflow-release.dat", trel)):
        if name == "sflow-release.dat" and c == "instant":
            continue
        with open(os.path.join(d, name), "w") as f:
            f.write("%d %d %.3f %.3f %.3f %.3f\n" % (len(x), len(y), x[0], y[0], dx, dx))
            for j in range(len(y)):
                f.write(" ".join("%.3f" % q for q in v[j]) + "\n")
print("block volume %.1f m3 (each half %.1f m3)" % (inblock.sum()*dx*dx*3.0, inblock.sum()*dx*dx*1.5))
