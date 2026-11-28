#!/usr/bin/env python3
# The fit of the agent body (tools/fish_body, plan 2 work package E), no REEF3D run needed:
#   * body2d.py with its extensions off is the body of fish_agents.py (body_step, n_up, n_down): the same speed
#   * fit_body.py recovers the constants of a body from beats made with them (fixed gaits from rest, a slow-down, a
#     coast, a brake, turns; the data with a noise of 0.005 m/s and 0.01 rad/s), with the relaxation in tail beats and
#     with the drag law (C_d)
#   * the drag law coasts as the analytic glide U0 / (1 + C_d U0 t / L)
#   * fish_agents.py with the defaults of body2d.py ("model": 2) makes the track of its own body
#   python3 fish_body_check.py        exit 1 on a failure
# Architect: Hans Bihs
import os, sys
import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..", "..", "tools", "fish_body"))
sys.path.insert(0, os.path.join(here, "..", "..", "tools", "fish_agents"))
import body2d as bd           # noqa: E402
import fit_body as fb         # noqa: E402
import fish_agents as fa      # noqa: E402

nfail = npass = 0


def check(ok, what):
    global nfail, npass
    print(("  ok    " if ok else "  FAIL  ") + what)
    if ok:
        npass += 1
    else:
        nfail += 1


# 1. body2d without its extensions is body_step of fish_agents.py
p = dict(bd.DEFAULT)
body = dict(fa.BODY, n_up=p["n_up"], n_down=p["n_down"], s=p["s"], tau_coast=p["tau_coast"], C_brake=p["C_brake"],
            k=p["k"], k_brake=p["k_brake"])
t = np.arange(0.0, 14.0, 0.01)
f = np.where(t < 6, 1.0, 0.5); a = np.where((t > 8) & (t < 10), 0.0, 1.0); b = np.where(t > 10, 2.5 * (t < 12), 0.0)
U, W, V = bd.simulate(t, f, a, b, p)
Ug, err = 0.0, 0.0
for i in range(1, len(t)):
    Ug, _ = fa.body_step(Ug, f[i - 1], a[i - 1], b[i - 1], t[i - 1], 0.01, body)
    err = max(err, abs(body["s"] * Ug - U[i]))
check(err < 1e-12, "body2d = body_step of fish_agents.py: largest difference of the speed %.1e m/s" % err)

# 2. the fit recovers the constants of a body from its beats
true = dict(bd.DEFAULT, s=0.55, s_f=0.12, s_b=0.3, n_up=0.9, n_down=0.7, tau_coast=9.0, C_brake=3.0, k=1.2, k_b=0.25,
            k_brake=1.6, beta_b=-0.15)
rng = np.random.default_rng(3)


def synth(name, T, g):
    t = np.arange(0.0, T, 0.01)
    f, a, b = (np.array([g(x)[k] for x in t]) for k in range(3))
    U, Wm, V = bd.simulate(t, f, a, b, true)
    ph = np.cumsum(2 * np.pi * f * 0.01); k = np.floor(ph / (2 * np.pi))
    edges = [0] + list(np.where(np.diff(k) > 0)[0] + 1) + [len(t) - 1]
    beats = []
    for i0, i1 in zip(edges[:-1], edges[1:]):
        if i1 - i0 < 5:
            continue
        sl = slice(i0, i1)
        beats.append({"t0": t[i0], "t1": t[i1], "U": U[sl].mean() + rng.normal(0, 0.005), "W": -Wm[sl].mean() + rng.normal(0, 0.01),
                      "V": V[sl].mean() + rng.normal(0, 0.005), "f": f[sl].mean(), "a": a[sl].mean(), "P": 0.0, "use": True})
    return {"name": name, "t": t, "f": f, "a": a, "b": b, "beats": beats}


def recover(true, free, label):
    rs = [synth(*r) for r in RUNS]          # made with the global true
    best, _ = fb.fit(free, dict(bd.DEFAULT, C_d=0.2 if "C_d" in free else 0.0), rs)
    # (k and k_b, and s and s_b, trade within the bends of the runs: compared as the turn constant and the stride there)
    rel = {k: abs(best[k] - true[k]) / max(abs(true[k]), 0.1) for k in free if k not in ("k_b", "s_b")}
    for bb in (0.3, 0.6):
        rel["k_eff(b %.1f)" % bb] = abs(best["k"] * (1 - best["k_b"] * bb) / (true["k"] * (1 - true["k_b"] * bb)) - 1.0)
        rel["s_eff(b %.1f)" % bb] = abs(bd.stride(dict(true, **best), 1.0, bb) / bd.stride(true, 1.0, bb) - 1.0)
    worst = max(rel, key=rel.get)
    check(rel[worst] < 0.1, "fit_body recovers the constants from noisy beats, %s (k, k_b and s, s_b as the turn constant "
          "and the stride in the bends): largest error %.1f %% (%s)" % (label, 100 * rel[worst], worst))


RUNS = [("f05", 14, lambda x: (0.5, 1, 0)), ("f075", 12, lambda x: (0.75, 1, 0)), ("f10", 10, lambda x: (1.0, 1, 0)),
        ("f15", 7, lambda x: (1.5, 1, 0)), ("down", 14, lambda x: (1.0 if x < 6 else 0.5, 1, 0)),
        ("coast", 14, lambda x: (1.0, 1.0 if x < 6 else 0.0, 0)), ("brake", 10, lambda x: (1.0, 1.0 if x < 6 else 0.0, 0 if x < 6 else 2.5)),
        ("t03", 8, lambda x: (1.0, 1, 0.3 if x > 4 else 0)), ("t06", 8, lambda x: (1.0, 1, 0.6 if x > 4 else 0)),
        ("t05", 11, lambda x: (0.5, 1, 0.6 if x > 6 else 0))]
recover(true, [k for k in bd.PARAMS if k != "C_d"], "relaxation in tail beats")
true = dict(true, C_d=0.45)
recover(true, [k for k in bd.PARAMS if k != "tau_coast"], "drag law of the glide")

# 2b. the drag law coasts as the analytic glide
pc = dict(bd.DEFAULT, C_d=0.4)
tc = np.arange(0.0, 8.0, 0.001)
Uc, _, _ = bd.simulate(tc, np.ones_like(tc), np.zeros_like(tc), np.zeros_like(tc), pc, U0=0.6)
ea = np.abs(Uc - 0.6 / (1 + 0.4 * 0.6 * tc)).max()
check(ea < 1e-3, "the drag law coasts as U0 / (1 + C_d U0 t / L): largest difference %.1e m/s" % ea)

# 3. fish_agents.py with "model": 2 and the defaults of body2d: the track of the old body (no brain decisions differ)
case = {"T": 8.0, "dt": 0.01, "out_dt": 0.5, "box": [0, 30, 0, 30], "walls": [0, 0, 0, 0], "flow": {"type": "uniform", "u": [0.0, 0.0]},
        "brain": {"U_des": 0.5}, "fish": [{"x": 2.0, "z": 15.0, "h": 0.0, "target": [20.0, 15.0]}]}
_, r1 = fa.run(dict(case, body=dict(n_up=0.8, n_down=0.65, s=0.5)))
_, r2 = fa.run(dict(case, body=dict(model=2, s=0.5)))
d = max(abs(x[2] - y[2]) + abs(x[3] - y[3]) for x, y in zip(r1, r2))
check(d < 1e-6, "fish_agents.py, the body of body2d.py (model 2) against its own body: largest difference of the track %.1e m" % d)

print("%d passed, %d failed" % (npass, nfail))
sys.exit(1 if nfail else 0)
