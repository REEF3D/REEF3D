#!/usr/bin/env python3
# Architect: Hans Bihs
"""
waverecon.dat (A [m], omega [rad/s], epsilon [rad]) for B 92 51 from the measured Mase & Kirby (1992)
record at the slope toe (refdata/waves/mase_kirby/mase_kirby_eta_h470mm.dat, eta in cm, dt = 0.05 s),
the same decomposition as FUNWAVE-TVD's fft4wavemaker.m: FFT of the first 8192 samples (409.6 s,
mean removed), components between 0.2 and 3.0 Hz. REEF3D's linear irregular waves use
eta = sum A cos(k x - omega t - epsilon) with x measured from the wave origin (B 105), so
epsilon = arg(X_n) reproduces the record at the origin: eta(0, t) = sum A cos(omega t + epsilon).

usage: mase_kirby_waverecon.py [fmin fmax]
"""
import cmath
import math
import os
import sys

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, ".."))
from benchmark import fft  # noqa: E402

fmin = float(sys.argv[1]) if len(sys.argv) > 1 else 0.2
fmax = float(sys.argv[2]) if len(sys.argv) > 2 else 3.0
fn = os.path.join(here, "..", "refdata", "waves", "mase_kirby", "mase_kirby_eta_h470mm.dat")
eta = [float(l.split()[0]) * 0.01 for l in open(fn) if l.strip() and l[0] != "#"]
dt, n = 0.05, 8192
x = eta[:n]
m = sum(x) / n
X = fft([complex(v - m) for v in x])
T = n * dt
with open("waverecon.dat", "w") as f:
    for k in range(1, n // 2):
        fr = k / T
        if fmin < fr < fmax:
            f.write("%.8e %.8f %.8f\n" % (2 * abs(X[k]) / n, 2 * math.pi * fr, cmath.phase(X[k])))
