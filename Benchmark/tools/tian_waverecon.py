#!/usr/bin/env python3
# Architect: Hans Bihs
"""
waverecon.dat (A [m], omega [rad/s], epsilon [rad]) for B 92 51: the constant-steepness deep-water focused wave packets
of Tian, Perlin & Choi (2010) as simulated by Derakhti & Kirby (2016, JFM 790, Table 1): depth 0.62 m, f_c = 1.7 Hz,
df/f_c = 0.824, N = 128 components uniformly spaced over f_c +- df/2, a_i = S/(k_i N), linear focus x_b = 7.0 m,
t_b = 25.0 s. REEF3D's linear irregular waves use eta = sum A cos(k x - omega t - epsilon), so
epsilon = k x_b - omega t_b focuses the components at (x_b, t_b).

usage: tian_waverecon.py S
"""
import math
import sys

S = float(sys.argv[1])
d, fc, dfr, N, xb, tb, g = 0.62, 1.7, 0.824, 128, 7.0, 25.0, 9.81


def kdisp(w):
    k = w * w / g
    for _ in range(200):
        k = w * w / (g * math.tanh(k * d))
    return k


df = dfr * fc
with open("waverecon.dat", "w") as f:
    for n in range(N):
        fn = fc - df / 2 + df * n / (N - 1)
        w = 2 * math.pi * fn
        k = kdisp(w)
        f.write("%.10e %.10e %.10e\n" % (S / (k * N), w, k * xb - w * tb))
