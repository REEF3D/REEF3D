# T2: lake at rest with an island (wetting and drying) on a wavy, skewed curvilinear grid:
# the velocities must stay at round-off
import sys, numpy as np
from p1lib import *

n, dx = 60, 0.5
control = """B 1 %g
B 10 0.0 %g 0.0 %g 0.0 0.5

M 10 1
M 20 2
""" % (dx, n*dx, n*dx)

ctrl = """A 10 2
A 210 3
A 211 4
A 220 0
A 240 1

F 60 1.0

N 41 100.0
N 47 0.3

M 10 1

W 22 -9.81
G 50 2
"""
steps = int(sys.argv[1]) if len(sys.argv) > 1 else 200

I, J = np.meshgrid(np.arange(n+1), np.arange(n+1), indexing="ij")
x, y = I*dx, J*dx
L = n*dx
# wavy, skewed nodes (boundary nodes move along the boundary only in the interior direction)
X = x + 0.25*y + 0.8*np.sin(2*np.pi*x/L)*np.sin(np.pi*y/L)
Y = y + 0.1*x + 0.6*np.sin(np.pi*x/L)*np.sin(4*np.pi*y/L)
rng = np.random.default_rng(1)
Z = 1.3*np.exp(-((X-17.0)**2+(Y-17.5)**2)/12.0) + 0.05*rng.standard_normal(X.shape)

for name, grid in [("wavy", (X, Y, Z))]:
    d = os.path.join(S, "runs/t2_"+name)
    rc = run(d, control, ctrl, grid=grid, steps=steps)
    st = final_state(d)
    wet = st["WL"] > 0.002
    print("%s rc=%d count=%d t=%.3f wet cells %d of %d" % (name, rc, st["count"], st["simtime"], wet.sum(), wet.size))
    print("  max |P| %.3e  max |Q| %.3e  max |eta| (wet) %.3e" % (np.abs(st["P"]).max(), np.abs(st["Q"]).max(), np.abs(st["eta"][wet]).max()))
    for l in log_lines(d, "geometry|curvilinear|!!!"):
        print("  ", l)
