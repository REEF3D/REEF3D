#!/usr/bin/env python3
# The constants of the agent body (fish_agents.py, body_step) fitted to a resolved 2D run of REEF3D.
#   python3 calibrate.py <REEF3D_6DOF_fish_0.dat> [t0 t1 [out.csv]] [--C_brake 3.5]
#     its own gait (f, a, b from the fish log) drives the body of fish_agents.py, and s, tau_up, tau_down and
#     k_brake are fitted to the speed and the heading of the resolved fish (least squares, Nelder-Mead); tau_coast
#     11 s and k 1.3 for a swimming turn (the 2D swimmer, Fish Swimming validations 13 and 08) and C_brake are
#     held; out.csv: t, U and heading of the run and the model
#   python3 calibrate.py <REEF3D_6DOF_fish_0.dat> --brake t0 t1
#     C_brake from a run that brakes from t0 (tail still, body bent): U = U0/(1 + C_brake U0 (t - t0))
#     (the 2D swimmer, b = 2.5 from 0.53 m/s, validation 16: 3.5 1/m, rms 0.018 m/s)
# Architect: Hans Bihs
import sys, math
import os
import numpy as np
from scipy.optimize import minimize
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fish_agents as fa                                        # noqa: E402


def table(p):
    rows = []
    for l in open(p):
        try:
            rows.append([float(x) for x in l.split()])
        except ValueError:
            pass
    n = max(set(len(r) for r in rows), key=lambda k: sum(len(r) == k for r in rows))
    return np.array([r for r in rows if len(r) == n])


args = [x for x in sys.argv[1:]]
CB = 3.5
if "--C_brake" in args:
    k = args.index("--C_brake"); CB = float(args[k + 1]); del args[k:k + 2]
f = table(args[0])
if "--brake" in args:
    k = args.index("--brake"); b0, b1 = float(args[k + 1]), float(args[k + 2])
    mb = (f[:, 0] >= b0) & (f[:, 0] <= b1)
    tb, Ub = f[mb, 0] - b0, f[mb, 1]
    fit = minimize(lambda q: np.mean((Ub[0] / (1.0 + q[0] * Ub[0] * tb) - Ub) ** 2), [3.0], method="Nelder-Mead")
    print("C_brake %.2f 1/m from %.3f m/s, rms %.4f m/s" % (fit.x[0], Ub[0], math.sqrt(fit.fun)))
    sys.exit(0)
sys.argv = [sys.argv[0]] + args
t0 = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
t1 = float(sys.argv[3]) if len(sys.argv) > 3 else f[-1, 0]
t, U, H = f[:, 0], f[:, 1], -np.radians(f[:, 3])          # 2D: the heading is minus the pitch angle
F, A, B = f[:, 9], f[:, 10], f[:, 11]
m = (t >= t0) & (t <= t1)


def model(p):
    body = dict(fa.BODY, s=p[0], tau_up=p[1], tau_down=p[2], C_brake=CB, k_brake=p[3], L=1.0)
    tt = np.arange(0.0, t1, 0.01)
    Ug, h = 0.0, H[0]
    Us, hs = [], []
    for x in tt:
        fg, ag, bg = np.interp(x, t, F), np.interp(x, t, A), np.interp(x, t, B)
        Ug, kt = fa.body_step(Ug, fg, ag, bg, x, 0.01, body)
        h += 0.01 * kt * bg * body["s"] * Ug
        Us.append(body["s"] * Ug); hs.append(h)
    return np.interp(t, tt, Us), np.interp(t, tt, hs)


def cost(p):
    if min(p) <= 0:
        return 1e9
    Um, hm = model(p)
    return np.mean((Um[m] - U[m]) ** 2) / 0.01 + np.mean((np.unwrap(hm[m]) - np.unwrap(H[m])) ** 2) / 0.25


p0 = [0.48, 0.6, 0.6, 1.8]
r = min((minimize(cost, x0, method="Nelder-Mead", options={"maxiter": 3000, "xatol": 1e-4, "fatol": 1e-7})
         for x0 in (p0, [0.5, 0.8, 0.4, 1.5])), key=lambda q: q.fun)
Um, hm = model(r.x)
print("s %.3f tau_up %.2f tau_down %.2f k_brake %.2f (C_brake %.1f 1/m held) | rms U %.3f m/s, rms heading %.1f deg" %
      (*r.x, CB, math.sqrt(np.mean((Um[m] - U[m]) ** 2)), math.degrees(math.sqrt(np.mean((np.unwrap(hm[m]) - np.unwrap(H[m])) ** 2)))))
if len(sys.argv) > 4:
    np.savetxt(sys.argv[4], np.c_[t[m], U[m], Um[m], np.degrees(H[m]), np.degrees(hm[m])], fmt="%.5f",
               header="t U_cfd U_model heading_cfd_deg heading_model_deg", comments="")
