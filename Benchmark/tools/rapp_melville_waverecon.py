#!/usr/bin/env python3
# Architect: Hans Bihs
"""
waverecon.dat (A [m], omega [rad/s], epsilon [rad]) for B 92 51: the constant-amplitude focused wave packets of
Rapp & Melville (1990) as simulated by Derakhti & Kirby (2016, JFM 790, Table 1): depth 0.6 m, f_c = 0.88 Hz,
df/f_c = 0.73, N = 32 components uniformly spaced over f_c +- df/2, a = S/(k_c N), linear focus x_b = 8.5 m,
t_b = 20.5 s. REEF3D's linear irregular waves use eta = sum A cos(k x - omega t - epsilon), so
epsilon = k x_b - omega t_b focuses the components at (x_b, t_b).

usage: rapp_melville_waverecon.py S
"""
import math
import sys

S = float(sys.argv[1])
d, fc, dfr, N, xb, tb, g = 0.6, 0.88, 0.73, 32, 8.5, 20.5, 9.81


def kdisp(w):
    k = w * w / g
    for _ in range(200):
        k = w * w / (g * math.tanh(k * d))
    return k


kc = kdisp(2 * math.pi * fc)
a = S / (kc * N)
df = dfr * fc
with open("waverecon.dat", "w") as f:
    for n in range(N):
        fn = fc - df / 2 + df * n / (N - 1)
        w = 2 * math.pi * fn
        k = kdisp(w)
        f.write("%.10e %.10e %.10e\n" % (a, w, k * xb - w * tb))
