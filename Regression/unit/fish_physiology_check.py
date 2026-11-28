#!/usr/bin/env python3
# The physiology of the agent fish (tools/fish_physiology/physiology.py) and the cage stages of fish_agents3d.py, no
# REEF3D run needed:
#   * U_crit(W, T): the fit of Hvas et al. (2021) at the fish of Athammer et al. (2024), held beyond its maximum in W
#   * below the aerobic limit no debt; at rest the oxygen uptake is SMR
#   * at a constant speed above it the fish is exhausted after C/(m - 1); a debt is paid back as the model says
#   * the U_crit test on the model gives back U_crit (the calibration); the fluctuating currents of Athammer et al.:
#     at 80 and 100 % U_crit most fish endure 4 h, at 140 % none does
#   * fish_agents3d.py with "physiology" and "stages": in a weak current the school circles, in a current above every
#     fish's preferred speed every fish stands on it facing upstream, and exhausted fish end on the downstream net
#   * the glide by a drag law (C_glide); in 0.9 m/s the fish burst and glide, at a lower cost than steady swimming, and
#     hold their station; a fish pressed against the net escapes (patch 0044)
#   python3 fish_physiology_check.py        exit 1 on a failure
# Architect: Hans Bihs
import math, os, sys
import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..", "..", "tools", "fish_physiology"))
sys.path.insert(0, os.path.join(here, "..", "..", "tools", "fish_agents"))
import physiology as ph       # noqa: E402

nfail = npass = 0


def check(ok, what):
    global nfail, npass
    print(("  ok    " if ok else "  FAIL  ") + what)
    if ok:
        npass += 1
    else:
        nfail += 1


# 1. U_crit(W, T)
u = float(ph.ucrit_cm(837.0, 9.0))
check(abs(u - 90.9) < 0.2, "U_crit fit at 837 g, 9 C: %.1f cm/s (90.9 by hand; Athammer et al. measured 94.5)" % u)
check(abs(float(ph.ucrit_cm(5000.0, 12.0)) - float(ph.ucrit_cm(1640.0, 12.0))) < 0.05,
      "U_crit fit held at its maximum beyond 1640 g: %.1f cm/s at 5 kg" % float(ph.ucrit_cm(5000.0, 12.0)))

# 2. below the aerobic limit; at rest
P = ph.Physiology(3, 1.0, vary=False)
for _ in range(3600):
    st = P.step(np.array([0.0, 0.5, 0.97 * P.p["u_a"]]), 1.0)
check(np.all(P.D == 0.0) and abs(st["MO2"][0] - P.p["r_s"] * P.p["MMR"]) < 1e-9 and abs(st["m"][2] - 1.0) < 0.1,
      "below the aerobic limit no debt in 1 h; at rest MO2 = SMR (%.0f mg/kg/h); demand at 0.97 U_a %.2f MMR"
      % (st["MO2"][0], st["m"][2]))

# 3. exhausted at a constant speed after C/(m - 1); the recovery
P = ph.Physiology(1, 1.0, vary=False)
U = np.array([1.2])
m = float(P.demand(U)[0])
t = 0.0
while not P.exhausted[0]:
    P.step(U, 0.1); t += 0.1
te = P.C[0] / (m - 1.0)
check(abs(t - te) < 0.11, "exhausted at 1.2 U_crit after %.1f s (C/(m - 1) = %.1f s)" % (t, te))
check(np.all(P.u_max() == P.u_aerobic()), "the exhausted fish swims no faster than U_a (%.3f m/s)" % P.u_aerobic()[0])
D0 = P.D[0]
tr = P.p["tau_r"]
for _ in range(600):
    P.step(np.array([0.0]), 1.0)
# at rest: all the spare scope 1 - r_s while D > tau_r, then first order
sp = 1.0 - P.p["r_s"]
t1 = max(0.0, (D0 - tr) / sp)
De = (D0 - sp * 600.0) if 600.0 <= t1 else min(D0, tr) * math.exp(-sp * (600.0 - t1) / tr)
check(abs(P.D[0] - De) < 1e-3 * D0, "the debt after 10 min at rest: %.1f s at MMR (the model: %.1f)" % (P.D[0], De))
check(not P.exhausted[0] if P.D[0] <= P.p["f_resume"] * P.C[0] else P.exhausted[0],
      "exhausted until the debt is back to f_resume C (debt %.0f of %.0f)" % (P.D[0], P.C[0]))

# 4. the calibration: the U_crit test on the model, the fluctuating currents
uc, _ = ph.brett(n=300, ucrit=ph.ATHAMMER_UCRIT, dt=2.0, seed=3)
check(abs(uc.mean() / ph.ATHAMMER_UCRIT - 1.0) < 0.03, "the U_crit test on the model: %.3f m/s (the fish: %.3f)" % (uc.mean(), ph.ATHAMMER_UCRIT))
end = []
for k in (0, 1, 3):
    g = ph.ATHAMMER[k]
    te, _ = ph.fluctuating(g["peak"], g["min"], g["mean"], g["period"], n=300, dt=1.0, seed=3)
    end.append(float(np.isnan(te).mean()))
check(end[0] > 0.95 and end[1] > 0.7 and end[2] < 0.05,
      "fluctuating currents: endured 4 h at 80 %% %.2f, 100 %% %.2f, 140 %% %.2f (observed 1, 1, 0)" % tuple(end))

# 5. the cage stages
import fish_agents3d as f3     # noqa: E402
base = {"T": 120.0, "dt": 0.02, "out_dt": 1.0, "seed": 4,
        "body": {"L": 0.6, "f0": 1.6, "s": 0.55, "n_up": 3.0, "n_down": 4.3, "tau_coast": 9.0, "C_brake": 3.5, "r_fish": 0.1, "at_rest": 1},
        "brain": {"U_des": 0.25, "b_max": 1.0, "D": 1.5, "s0": 0.55, "n_u": 3.0, "n_d": 4.3, "s_sd": 0.1, "tau_c": 9.0},
        "box": [0, 20, 0, 16, 0, 7], "walls": [0, 0, 1, 1], "nets": [[10.0, 8.0, 4.0, 1.0, 6.0]],
        "physiology": {"W": 2500, "T": 12}, "stages": {"U_pref": [0.30, 0.03], "h": 0.1, "lead": 1.0},
        "fish": {"n": 20, "school": {"c": [10.0, 8.0], "r": [1.5, 3.3], "z": [2.0, 5.0], "U": 0.25, "dir": 1,
                 "d_sep": 0.5, "d_ali": 1.5, "d_coh": 2.0, "lead": 1.0, "leash": 1.8}}}
res = {}
for name, u, extra in (("weak", 0.05, {}), ("strong", 0.6, {}), ("exhaust", 0.95, {"C": 2.0, "cv_C": 0.0, "f_resume": 0.0})):
    c = dict(base, flow={"type": "uniform", "u": [u, 0.0, 0.0]}, physiology=dict(base["physiology"], **extra))
    fish, rows = f3.run(c)
    R = np.array(rows); R = R[R[:, 0] > 90.0]
    res[name] = R
Rw, Rs, Re = res["weak"], res["strong"], res["exhaust"]
check(np.all(Rw[:, 11] == 0), "cage, 0.05 m/s: every fish circles (stage 0 in %.0f %% of the samples)" % (100 * np.mean(Rw[:, 11] == 0)))
hs = np.radians(Rs[:, 5]); up = np.mean(np.cos(hs) < -0.9)
check(np.all(Rs[:, 11] >= 1) and up > 0.9, "cage, 0.6 m/s: every fish stands on the current (stage 1-2), %.0f %% facing upstream" % (100 * up))
xr = (Re[:, 2] - 10.0) / 4.0
check(np.all(Re[:, 11] == 3) and np.mean(xr > 0.7) > 0.9,
      "cage, 0.95 m/s, little capacity: every fish exhausted, %.0f %% on the downstream net (x > 0.7 r)" % (100 * np.mean(xr > 0.7)))


# 6. the burst and glide, the glide by a drag law, the escape from the net (patch 0044)
Ug, t, C, L_ = 1.0 / 0.55, 0.0, 0.13, 0.6
bd = dict(base["body"], C_glide=C, k=1.3, k_brake=1.77)
for _ in range(100):
    Ug, _k = f3.body_step3(Ug, 1.6, 0.0, 0.0, 1.0, 0.01, bd); t += 0.01
Ue = 1.0 / (1.0 + C * 1.0 * t / L_)
check(abs(0.55 * Ug - Ue) < 2e-3, "the glide by the drag law C_glide 0.13: %.4f m/s after 1 s (U0/(1 + C U0 t/L) = %.4f)" % (0.55 * Ug, Ue))
out = {}
for name, bg, cg in (("steady", 0, 0.0), ("bg", 1, 0.13)):
    c = dict(base, flow={"type": "uniform", "u": [0.9, 0.0, 0.0]}, body=dict(base["body"], C_glide=cg),
             stages=dict(base["stages"], burst_glide=bg))
    fish, rows = f3.run(c)
    R = np.array(rows); out[name] = R[R[:, 0] > 60.0]
Rst, Rbg = out["steady"], out["bg"]
gl = np.mean(Rbg[:, 8] < 0.5)
check(np.all(Rbg[:, 11] == 2) and 0.2 < gl < 0.8, "cage, 0.9 m/s: every fish bursts and glides (stage 2), gliding in %.0f %% of the samples" % (100 * gl))
check(Rbg[:, 13].mean() < 0.85 * Rst[:, 13].mean() and np.mean((Rbg[:, 2] - 10.0) / 4.0) < 0.4,
      "cage, 0.9 m/s: the burst and glide costs less (MO2 %.0f against %.0f steady) and holds the station (mean x/r %.2f)"
      % (Rbg[:, 13].mean(), Rst[:, 13].mean(), np.mean((Rbg[:, 2] - 10.0) / 4.0)))
f = f3.Fish3(0, 13.85, 8.0, 3.0, 0.0, None, dict(base["body"]), base["brain"], 0.0)     # on the downstream net, facing downstream
f.stage, f.station, f.goal, f.pin = 1, [11.0, 8.0, 3.0], [10.4, 8.0, 3.0], 2.5
st = {"net_escape": 1, "t_pin": 2.0, "net": [10.0, 8.0, 4.0, 1.0, 6.0]}
f3.gait_override(f, 0, 100.0, st, (np.array([0.5]), np.array([0.0])), dict(base["body"]), 0.6)
check(f.escape is not None and f.g["f"][2] == 3.2 and f.g["a"][2] == 1.0 and abs(f.g["b"][2]) == 1.0,
      "the escape: a fish pressed against the net for 2.5 s turns hard (b %.1f) at its fastest gait (%.1f Hz)" % (f.g["b"][2], f.g["f"][2]))

print("%d passed, %d failed" % (npass, nfail))
sys.exit(1 if nfail else 0)
