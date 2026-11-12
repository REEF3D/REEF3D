#!/usr/bin/env python3
# The constants of the agent body (fish_agents.py) fitted to a resolved 2D run of REEF3D: its own gait (f, a, b
# from REEF3D_6DOF_fish_<n>.dat) drives the body model of fish_agents.py, and s, tau_up, tau_down, tau_brake,
# U_brake, k_brake are fitted to the speed and the heading of the resolved fish (least squares, Nelder-Mead);
# tau_coast 11 s and k 1.3 for a swimming turn (the 2D swimmer, Fish Swimming validations 13 and 08) are held.
# The 2D swimmer at its station with the brain (validation 17, 0-12 s): s 0.50, tau_up 0.61 s, tau_down 0.49 s,
# tau_brake 0.87 s, U_brake 0.13 m/s, k_brake 1.77; rms 0.04 m/s and 5 deg.
#   python3 calibrate.py <REEF3D_6DOF_fish_0.dat> [t0 t1 [out.csv]]     out.csv: t, U and heading of the run and the model
# Architect: Hans Bihs
import sys, math
import numpy as np
from scipy.optimize import minimize


def table(p):
    rows = []
    for l in open(p):
        try:
            rows.append([float(x) for x in l.split()])
        except ValueError:
            pass
    n = max(set(len(r) for r in rows), key=lambda k: sum(len(r) == k for r in rows))
    return np.array([r for r in rows if len(r) == n])


f = table(sys.argv[1])
t0 = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
t1 = float(sys.argv[3]) if len(sys.argv) > 3 else f[-1, 0]
t, U, H = f[:, 0], f[:, 1], -np.radians(f[:, 3])          # 2D: the heading is minus the pitch angle
F, A, B = f[:, 9], f[:, 10], f[:, 11]
m = (t >= t0) & (t <= t1)


def model(p):
    s, tau_up, tau_down, tau_brake, U_brake, k_brake = p
    tau_coast, k = 11.0, 1.3                           # the coasting time (validation 13) and k of a swimming turn (validation 08)
    tt = np.arange(0.0, t1, 0.01)
    Ug, h = 0.0, H[0]
    Us, hs = [], []
    for x in tt:
        fg, ag, bg = np.interp(x, t, F), np.interp(x, t, A), np.interp(x, t, B)
        target = fg * ag * min(1.0, x / 1.0)
        if ag < 0.5:
            if abs(bg) > 0.5:
                tau = tau_brake; target = max(target, U_brake / s)
            else:
                tau = tau_coast
        else:
            tau = tau_up if target >= Ug else tau_down
        Ug += 0.01 / tau * (target - Ug)
        h += 0.01 * (k_brake if ag < 0.5 else k) * bg * s * Ug
        Us.append(s * Ug); hs.append(h)
    return np.interp(t, tt, Us), np.interp(t, tt, hs)


def cost(p):
    if min(p) <= 0:
        return 1e9
    Um, hm = model(p)
    return np.mean((Um[m] - U[m]) ** 2) / 0.01 + np.mean((np.unwrap(hm[m]) - np.unwrap(H[m])) ** 2) / 0.25


p0 = [0.48, 0.8, 0.8, 1.9, 0.12, 1.3]
r = minimize(cost, p0, method="Nelder-Mead", options={"maxiter": 1500, "xatol": 1e-3, "fatol": 1e-5})
Um, hm = model(r.x)
print("s %.3f tau_up %.2f tau_down %.2f tau_brake %.2f U_brake %.3f k_brake %.2f | rms U %.3f m/s, rms heading %.1f deg" %
      (*r.x, math.sqrt(np.mean((Um[m] - U[m]) ** 2)), math.degrees(math.sqrt(np.mean((np.unwrap(hm[m]) - np.unwrap(H[m])) ** 2)))))
if len(sys.argv) > 4:
    np.savetxt(sys.argv[4], np.c_[t[m], U[m], Um[m], np.degrees(H[m]), np.degrees(hm[m])], fmt="%.5f",
               header="t U_cfd U_model heading_cfd_deg heading_model_deg", comments="")
