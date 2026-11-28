#!/usr/bin/env python3
# The fit of the fatigue model of physiology.py to the fluctuating currents of Athammer et al. (2024):
#   the U_crit test of the same fish (+0.15 m/s every 30 min from 0.15 m/s) must give back their U_crit (0.945 m/s);
#   80 and 100 % U_crit peaks (1 min): all endure 4 h; 120 %: 4 of 18 endure, the others are exhausted after 119 min
#   on average; 140 %: all exhausted, after 51 min on average; and a steady 80 % U_crit is sustained for 4 h (Hvas &
#   Oppedal 2017).
# The 0.5 and 2 min periods at 120 % are left out of the fit: they check it (113 and 123 min, 3 and 2 of 18).
# Fitted: u_a (the aerobic limit / U_crit), C (the anaerobic capacity, s at MMR), k_r = tau_r/C (the recovery), cv_C;
# with --free also cv_ucrit (the spread of U_crit, else fixed at --cv_ucrit) and gamma (the power-law demand).
#   python3 calibrate.py [--n 300] [--starts 1] [--free cv_ucrit] [--out fit.json]
# The constants of physiology.py come from "--n 300 --starts 1 --free cv_ucrit" (validation 32, about 15 min).
# Architect: Hans Bihs
import argparse, json, math, os, sys
import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import physiology as ph        # noqa: E402

FIT = [0, 1, 2, 3]             # the groups of ph.ATHAMMER in the fit


MODEL = {}                     # fixed constants (--gamma, --cv_ucrit) and the extra free ones (--free)


def unpack(q):
    kw = dict(MODEL.get("fixed", {}), u_a=q[0], C=math.exp(q[1]), tau_r=math.exp(q[1]) * math.exp(q[2]), cv_C=abs(q[3]))
    for j, k in enumerate(MODEL.get("free", [])):
        kw[k] = abs(q[4 + j])
    return kw


def evaluate(kw, n=400, groups=(0, 1, 2, 3, 4, 5, "steady"), dt=1.0, seed=7):
    uc, _ = ph.brett(n=n, ucrit=ph.ATHAMMER_UCRIT, dt=2.0, seed=seed, **kw)
    res = {"brett": float(uc.mean()), "brett_sd": float(uc.std()), "groups": []}
    if "steady" in groups or groups == FIT:
        te, _ = ph.steady(ph.SUSTAINED["u_rel"] * ph.ATHAMMER_UCRIT, ph.SUSTAINED["T"], ucrit=ph.ATHAMMER_UCRIT, n=n, dt=2.0 * dt, seed=seed, **kw)
        res["steady"] = {"name": ph.SUSTAINED["name"], "endured": float(np.isnan(te).mean()), "endured_obs": 1.0}
    for k in [g for g in groups if g != "steady"]:
        g = ph.ATHAMMER[k]
        te, _ = ph.fluctuating(g["peak"], g["min"], g["mean"], g["period"], ucrit=ph.ATHAMMER_UCRIT, n=n, dt=dt, seed=seed, **kw)
        fat = ~np.isnan(te)
        res["groups"].append({"name": g["name"], "endured": float(1 - fat.mean()), "t_fat": float(te[fat].mean() / 60) if fat.any() else None,
                              "t_fat_sd": float(te[fat].std() / 60) if fat.sum() > 1 else None,
                              "endured_obs": g["endured"], "t_fat_obs": g["t_fat"]})
    return res


def cost(q, n):
    if not (0.5 < q[0] < 1.05) or abs(q[3]) > 1.5 or not (-4.0 < q[2] < 4.0):
        return 1e3
    for j, k in enumerate(MODEL.get("free", [])):
        if not (MODEL["bounds"][k][0] <= abs(q[4 + j]) <= MODEL["bounds"][k][1]):
            return 1e3
    r = evaluate(unpack(q), n=n, groups=FIT, dt=2.0)
    J = 25.0 * (r["brett"] / ph.ATHAMMER_UCRIT - 1.0) ** 2 * 100
    J += 4.0 * (r["steady"]["endured"] - 1.0) ** 2
    for g in r["groups"]:
        J += 4.0 * (g["endured"] - g["endured_obs"]) ** 2
        if g["t_fat_obs"] is not None:
            tf = g["t_fat"] if g["t_fat"] is not None else 240.0
            J += ((tf - g["t_fat_obs"]) / 30.0) ** 2
    return J


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=400)
    ap.add_argument("--out", default="fit_physiology.json")
    ap.add_argument("--gamma", type=float, default=0.0, help="the demand as SMR plus a power gamma of the speed (0: exponential)")
    ap.add_argument("--cv_ucrit", type=float, default=0.07)
    ap.add_argument("--free", default="", help="also fit: gamma, cv_ucrit (comma separated)")
    ap.add_argument("--starts", type=int, default=3, help="the starts of the simplex (1-3)")
    a = ap.parse_args()
    MODEL["fixed"] = {"gamma": a.gamma, "cv_ucrit": a.cv_ucrit}
    MODEL["free"] = [k for k in a.free.split(",") if k]
    MODEL["bounds"] = {"gamma": (1.0, 6.0), "cv_ucrit": (0.0, 0.10)}
    start = {"gamma": 3.0, "cv_ucrit": 0.04}
    best = None
    for q0 in ([0.80, math.log(1200), math.log(0.4), 0.05], [0.80, math.log(1500), math.log(0.2), 0.5],
               [0.90, math.log(600), math.log(1.0), 0.2])[:a.starts]:
        q0 = q0 + [start[k] for k in MODEL["free"]]
        nq = len(q0)
        r = minimize(cost, q0, args=(a.n,), method="Nelder-Mead", options={"maxfev": 300, "xatol": 1e-3, "fatol": 1e-3,
                     "initial_simplex": np.array(q0) + np.vstack([np.zeros(nq), np.diag([0.04, 0.3, 0.3, 0.15] + [0.3 * start[k] for k in MODEL["free"]])])})
        print("start %s: J %.4f at %s" % (np.round(q0, 3), r.fun, {k: round(v, 3) for k, v in unpack(r.x).items()}), flush=True)
        if best is None or r.fun < best.fun:
            best = r
    kw = unpack(best.x)
    res = evaluate(kw, n=2000, dt=1.0)
    res["params"] = kw; res["J"] = best.fun
    json.dump(res, open(a.out, "w"), indent=1)
    print(json.dumps(res, indent=1))
