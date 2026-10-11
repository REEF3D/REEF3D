# Validation case cavitation 01: collapse of a vapour layer of thickness h0 at the closed end of a
# channel of length L under a liquid column, outlet pressure p_out, vapour pressure p_v.
# The liquid column between the layer and the outlet moves as a piston:
#
#     rho (L - h) d2h/dt2 = -(p_out - p_v),   h(0) = h0, dh/dt(0) = 0
#
#   python slab_compare.py [REEF3D-CFD-Cavitation.dat] [--plot]     in the run directory
# Reads rho (W 1), p_v (W 61), the layer (W 67) and the channel (B 10 in control.txt) from the inputs;
# the layer thickness of the CFD is the vapour volume over the channel cross-section.
import sys, os
import numpy as np

fn = next((a for a in sys.argv[1:] if not a.startswith("--")), "REEF3D_CFD_Cavitation/REEF3D-CFD-Cavitation.dat")
rho, pv, pout = 1000.0, None, 0.0
for l in open("ctrl.txt"):
    s = l.split()
    if len(s) >= 3 and s[0] == "W" and s[1] == "1": rho = float(s[2])
    if len(s) >= 3 and s[0] == "W" and s[1] == "61": pv = float(s[2])
    if len(s) >= 8 and s[0] == "W" and s[1] == "67": h0 = float(s[3]) - float(s[2])
for l in open("control.txt"):
    s = l.split()
    if len(s) == 8 and s[0] == "B" and s[1] == "10":
        x0, x1, y0, y1, z0, z1 = map(float, s[2:])
L, A = x1 - x0, (y1 - y0)*(z1 - z0)

d = np.loadtxt(fn, skiprows=1)
t, h = d[:, 0], d[:, 1]/A

def f(y):
    return np.array([y[1], -(pout - pv)/(rho*(L - y[0]))])
y, tt, dt, T, H = np.array([h0, 0.0]), 0.0, 1.0e-5, [0.0], [h0]
while y[0] > 0.0:
    k1 = f(y); k2 = f(y + dt/2*k1); k3 = f(y + dt/2*k2); k4 = f(y + dt*k3)
    y = y + dt/6*(k1 + 2*k2 + 2*k3 + k4); tt += dt
    T.append(tt); H.append(max(y[0], 0.0))
T, H = np.array(T), np.array(H)

tc = t[np.argmax(h < 1.0e-3*h0)] if np.any(h < 1.0e-3*h0) else float("nan")
print(f"L = {L} m, h0 = {h0} m, p_out - p_v = {pout - pv} Pa, rho = {rho}")
print(f"collapse time: piston {T[-1]:.4f} s, CFD {tc:.4f} s  ({100*(tc/T[-1] - 1):+.1f} %)")
print("    t [s]   h piston [m]   h CFD [m]")
for q in np.linspace(0.0, T[-1], 9):
    print(f"  {q:7.4f}   {np.interp(q, T, H):10.4f}   {np.interp(q, t, h):10.4f}")
if "--plot" in sys.argv:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(T, H, "k", lw=1.5, label="piston equation", zorder=3)
    ax.plot(t[::10], h[::10], "o", ms=4, mfc="none", label="REEF3D, Schnerr-Sauer")
    ax.set_xlabel("t [s]")
    ax.set_ylabel("vapour layer thickness h [m]")
    ax.legend()
    fig.tight_layout()
    fig.savefig("slab_collapse.png", dpi=130)
    print("written", os.path.abspath("slab_collapse.png"))
