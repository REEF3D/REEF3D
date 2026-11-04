# T3: 90 deg bend (centreline radius 50 m, width 20 m) with straight inlet and outlet reaches,
# flat bed, 2 m depth, 20 m3/s; steady state: discharge through the outlet, superelevation in the bend
import sys, numpy as np
from p1lib import *

R, B, Lin, Lout, ds, nj = 50.0, 20.0, 40.0, 60.0, 1.0, 20
Lb = 0.5*np.pi*R
ni = int(round((Lin+Lb+Lout)/ds))
s = np.linspace(-Lin, Lb+Lout, ni+1)
nrm = np.linspace(-0.5*B, 0.5*B, nj+1)          # from the right bank (outer bank of a left bend)

def centre(sv):
    if sv <= 0.0: return np.array([sv, 0.0]), np.array([0.0, 1.0])
    if sv <= Lb:
        t = sv/R
        return np.array([R*np.sin(t), R-R*np.cos(t)]), np.array([-np.sin(t), np.cos(t)])
    return np.array([R, R+(sv-Lb)]), np.array([-1.0, 0.0])

X = np.zeros((ni+1, nj+1)); Y = np.zeros_like(X)
for i, sv in enumerate(s):
    c, nl = centre(sv)
    X[i, :] = c[0] + nrm*nl[0]
    Y[i, :] = c[1] + nrm*nl[1]
Z = np.zeros_like(X)

np_ = int(sys.argv[2]) if len(sys.argv) > 2 else 1
control = """C 11 1
C 12 3
C 13 3
C 14 2

B 1 1.0
B 10 0.0 %g 0.0 %g 0.0 1.0

M 10 %d
M 20 2
""" % (ni, nj, np_)

ctrl = """A 10 2
A 210 3
A 211 4
A 220 0
A 240 1
A 218 1

B 50 0.05

B 60 1

F 50 2
F 60 2.0
F 62 2.0

N 41 400.0
N 47 0.3

M 10 %d

P 10 1
P 30 50.0

W 10 20.0
W 22 -9.81
G 50 2
""" % np_
T = float(sys.argv[1]) if len(sys.argv) > 1 else 400.0
ctrl = ctrl.replace("N 41 400.0", "N 41 %g" % T)

d = os.path.join(S, "runs/t3_bend_np%d" % np_)
rc = run(d, control, ctrl, grid=(X, Y, Z), np_=np_)
print("rc", rc)
q = log_lines(d, r"^(Inflow|Outflow):")
print(q[-2]); print(q[-1])
st = final_state(d)
eta = st["eta"].reshape(ni, nj)
imid = int(round((Lin+0.5*Lb)/ds))
print("bend middle (i=%d): eta right/outer bank %.5f, left/inner %.5f, difference %.5f" % (imid, eta[imid, 0], eta[imid, -1], eta[imid, 0]-eta[imid, -1]))
U = 20.0/(B*2.0)
print("estimate U^2 B/(g R) = %.5f" % (U*U*B/(9.81*R)))
print("free vortex: eta(n) across the bend middle:", " ".join("%.4f" % v for v in eta[imid, ::4]))
