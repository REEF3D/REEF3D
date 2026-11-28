#!/usr/bin/env python3
# Plan 2, work package E: the body of the agent fish (body2d.py) fitted to a set of resolved runs of REEF3D at once,
# with the uncertainty of each constant, and written as one body file.
#
#   python3 fit_body.py runs.json [--out body.json] [--ctrl body_ctrl.txt] [--boot 30] [--boot-mode beats|runs]
#                     [--free s,s_f,...] [--table fit.csv]
#
# runs.json: {"body": {"L": 1.0, "f0": 1.0, "t_ramp": 1.0}, "box": [x0, x1, z0, z1], "margin": 0.8,
#             "runs": [{"name": "f050", "dir": "<folder with the REEF3D_6DOF_* logs>", "current": [cx, cz]}, ...],
#             (a run may have its own "box"), "plane": "xz" (the 2D swimmer, default) | "xy" (a 3D body swimming in
#             the horizontal plane, its heading the yaw, e_x = (cos h, sin h); the box and the current in x, y)}
# Each run: its own gait (f, a, b of the fish log, every step) drives the model; per tail beat (the beats log of patch
# 0038) the mean speed through the water along the axis and the turn rate of the resolved body are the data (beats with
# the body within "margin" of the box are left out; the turn rate not before the end of the ramp of the gait), and the
# mean velocity across the axis through the water when the slip beta_b is among the free constants (by default it is
# not). By default the free constants are s, s_f, s_b, n_up, n_down, C_d, C_brake, k, k_b, k_brake: the body glides
# by the drag law (body2d.py; validation 30: the relaxation in tail beats cannot follow both the slower gait and the
# glide, and the glide is not exponential). The
# constants are fitted by least squares on all beats (the residuals in units of 0.01 m/s, 0.02 rad/s, 0.01 m/s); their
# uncertainty is the spread of the fits of resampled data (--boot n): by default the beats of each run are resampled
# (every run stays in, so the constants that one run alone sets, the coast and the brake, keep their data); with
# --boot-mode runs whole runs are resampled (the spread over gaits; a constant set by one run is then free in about a
# third of the samples). The output has the standard deviation and the 10-90 % range of each. The power per
# beat P = c3 (f a)^3 + c1 (f a) U^2 is fitted apart (linear). The heading of REEF3D turns with -omega of body2d.
#
# The body file: the constants, their uncertainty, the data and the fit per run; --ctrl writes the lines of the
# brain for REEF3D (X 200 s k, X 202 n_up n_down, X 203, X 204 tau_coast; with the drag law X 209 C_d and X 204 its
# exponential equivalent) from it, and fish_agents.py reads it
# ("body_file").
#
# Architect: Hans Bihs
import argparse, json, math, os, sys
import numpy as np
from scipy.optimize import least_squares

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import body2d as bd                                                   # noqa: E402

SIG = (0.01, 0.02, 0.01)       # the units of the residuals: speed [m/s], turn rate [rad/s], slip [m/s]


def table(p):
    rows = []
    for l in open(p):
        try:
            rows.append([float(x) for x in l.split()])
        except ValueError:
            pass
    n = max(set(len(r) for r in rows), key=lambda k: sum(len(r) == k for r in rows))
    return np.array([r for r in rows if len(r) == n])


def load(run, box, margin, plane="xz"):
    """the gait every step and the data per beat of a resolved run"""
    d = run["dir"]
    F = np.loadtxt(os.path.join(d, "REEF3D_6DOF_fish_0.dat"), skiprows=1)
    _, iu = np.unique(F[:, 0], return_index=True); F = F[iu]
    P = np.loadtxt(os.path.join(d, "REEF3D_6DOF_position_0.dat"), skiprows=1)
    _, ip = np.unique(P[:, 0], return_index=True); P = P[ip]
    Bt = table(os.path.join(d, "REEF3D_6DOF_fish_beats_0.dat"))
    t, Us, Ul, h = F[:, 0], F[:, 1], F[:, 2], np.unwrap(np.radians(F[:, 3]))
    cx, cz = run.get("current", [0.0, 0.0])
    # through the water, in the body frame (2D: e_x = (cos h, -sin h), e_l = (sin h, cos h) in x-z; "xy": e_x =
    # (cos h, sin h), e_l = (-sin h, cos h) in x-y)
    sg = -1.0 if plane == "xz" else 1.0
    Uw = Us - (cx * np.cos(h) + sg * cz * np.sin(h))
    Vw = Ul - (-sg * cx * np.sin(h) + cz * np.cos(h))
    kz = 3 if plane == "xz" else 2
    cum = lambda y: np.concatenate([[0.0], np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(t))])
    CU, CV = cum(Uw), cum(Vw)
    beats = []
    for r in Bt:
        t0, t1 = r[1], r[2]
        if t1 > t[-1] or t1 - t0 <= 0.0:
            continue
        x0, x1 = np.interp([t0, t1], P[:, 0], P[:, 1]); z0, z1 = np.interp([t0, t1], P[:, 0], P[:, kz])
        inside = min(x0, x1) > box[0] + margin and max(x0, x1) < box[1] - margin and \
            min(z0, z1) > box[2] + margin and max(z0, z1) < box[3] - margin
        T = t1 - t0
        beats.append({"t0": t0, "t1": t1, "U": (np.interp(t1, t, CU) - np.interp(t0, t, CU)) / T,
                      "W": (np.interp(t1, t, h) - np.interp(t0, t, h)) / T,
                      "V": (np.interp(t1, t, CV) - np.interp(t0, t, CV)) / T,
                      "f": r[3], "a": r[4], "P": r[8], "use": bool(inside)})
    if F.shape[1] >= 12:
        f, a, b = F[:, 9], F[:, 10], F[:, 11]
    else:
        # no gait control (X 199): the frequency of X 195 (of the run, else of its beats), a = 1, b = 0
        f0 = run.get("f", float(np.median(Bt[:, 3])))
        f, a, b = np.full_like(t, f0), np.ones_like(t), np.zeros_like(t)
    # the body at rest over the ground at t = 0 (as REEF3D starts it): in a current the water passes it, U0 = -c.e_x
    U0 = float(Uw[0])
    return {"name": run["name"], "t": t, "f": f, "a": a, "b": b, "beats": beats, "U0": U0}


def model_beats(R, p):
    """the model along the gait of the run, per beat: U, the turn rate (REEF3D's sense), the slip"""
    t = R["t"]
    U, W, V = bd.simulate(t, R["f"], R["a"], R["b"], p, R.get("U0", 0.0))
    cum = lambda y: np.concatenate([[0.0], np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(t))])
    CU, CW, CV = cum(U), cum(W), cum(V)
    out = []
    for b in R["beats"]:
        T = b["t1"] - b["t0"]
        out.append(((np.interp(b["t1"], t, CU) - np.interp(b["t0"], t, CU)) / T,
                    -(np.interp(b["t1"], t, CW) - np.interp(b["t0"], t, CW)) / T,
                    (np.interp(b["t1"], t, CV) - np.interp(b["t0"], t, CV)) / T))
    return out


def residuals(q, free, base, runs):
    p = dict(base); p.update(zip(free, q))
    res = []
    for R in runs:
        for b, m in zip(R["beats"], model_beats(R, p)):
            if b["use"]:
                res.append((m[0] - b["U"]) / SIG[0])
                # the turn rate not in the first beats (the yaw of the start, before the end of the ramp of the gait)
                if b["t0"] >= p["t_ramp"]:
                    res.append((m[1] - b["W"]) / SIG[1])
                # the slip only when it is fitted (the 2D swimmer crabs at 6-7 deg to its mean axis whatever its bias,
                # with the sign of its gait: no constant of the body, validation 30)
                if "beta_b" in free:
                    res.append((m[2] - b["V"]) / SIG[2])
    return np.array(res)


BOUNDS = {"s": (0.1, 1.5), "s_f": (-1.0, 1.0), "s_b": (-1.0, 1.0), "n_up": (0.1, 10.0), "n_down": (0.1, 20.0),
          "tau_coast": (0.5, 60.0), "C_d": (0.01, 10.0), "C_brake": (0.1, 30.0), "k": (-5.0, 5.0), "k_b": (-1.0, 1.0), "k_brake": (-5.0, 5.0),
          "beta_b": (-1.0, 1.0)}


def fit(free, base, runs):
    lo = [BOUNDS[k][0] for k in free]; hi = [BOUNDS[k][1] for k in free]
    q0 = [min(max(base[k], l + 1e-9 * (h - l)), h - 1e-9 * (h - l)) for k, l, h in zip(free, lo, hi)]
    r = least_squares(residuals, q0, bounds=(lo, hi), args=(free, base, runs), x_scale="jac")
    return dict(zip(free, r.x)), r


def fit_power(runs):
    rows = [(b["f"] * b["a"], b["U"], b["P"]) for R in runs for b in R["beats"] if b["use"] and b["a"] > 0.5]
    A = np.array([[(fa_) ** 3, fa_ * U * U] for fa_, U, _ in rows]); y = np.array([P for _, _, P in rows])
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    rms = float(np.sqrt(np.mean((A @ c - y) ** 2))) if len(y) else 0.0
    return {"c3": float(c[0]), "c1": float(c[1])}, rms, len(y)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs"); ap.add_argument("--out", default="body.json"); ap.add_argument("--ctrl", default="")
    ap.add_argument("--boot", type=int, default=30); ap.add_argument("--free", default="s,s_f,s_b,n_up,n_down,C_d,C_brake,k,k_b,k_brake")
    ap.add_argument("--table", default=""); ap.add_argument("--boot-mode", default="beats", choices=("beats", "runs"))
    a = ap.parse_args()
    cfg = json.load(open(a.runs))
    base = dict(bd.DEFAULT, **cfg.get("body", {}))
    runs = [load(r, r.get("box", cfg.get("box", [0, 8, 0, 6])), cfg.get("margin", 0.8), cfg.get("plane", "xz")) for r in cfg["runs"]]
    free = [k for k in a.free.split(",") if k]
    if "C_d" in free and base.get("C_d", 0.0) <= 0.0:
        base["C_d"] = 0.4                  # the start of the drag law (off by default in body2d.py)
    best, r = fit(free, base, runs)
    p = dict(base); p.update(best)
    # the bootstrap: the beats within each run, or whole runs
    rng = np.random.default_rng(1)
    boot = []
    for i in range(a.boot):
        if a.boot_mode == "runs":
            pick = [runs[j] for j in rng.integers(0, len(runs), len(runs))]
        else:
            pick = [dict(R, beats=[R["beats"][j] for j in rng.integers(0, len(R["beats"]), len(R["beats"]))]) for R in runs]
        try:
            bb, _ = fit(free, p, pick)
            boot.append([bb[k] for k in free])
        except Exception:
            pass
    boot = np.array(boot) if boot else np.zeros((0, len(free)))
    sd = {k: float(boot[:, i].std(ddof=1)) if len(boot) > 1 else None for i, k in enumerate(free)}
    ci = {k: [float(v) for v in np.percentile(boot[:, i], [10, 90])] if len(boot) > 1 else None for i, k in enumerate(free)}
    pw, prms, npw = fit_power(runs)
    p.update(pw)
    # the fit per run
    per = []
    for R in runs:
        m = model_beats(R, p)
        e = np.array([[mm[0] - b["U"], mm[1] - b["W"], mm[2] - b["V"]] for b, mm in zip(R["beats"], m) if b["use"]])
        per.append({"run": R["name"], "beats": int(len(e)),
                    "rms_U": float(np.sqrt((e[:, 0] ** 2).mean())) if len(e) else None,
                    "rms_W": float(np.sqrt((e[:, 1] ** 2).mean())) if len(e) else None,
                    "rms_V": float(np.sqrt((e[:, 2] ** 2).mean())) if len(e) else None})
    out = {"body": {k: (float(v) if isinstance(v, (int, float, np.floating)) else v) for k, v in p.items()},
           "sd": sd, "ci_10_90": ci, "boot_mode": a.boot_mode, "free": free, "boot": boot.tolist(), "power_rms": prms, "power_beats": npw, "runs": per,
           "turn_sign": "REEF3D heading rate = -omega of body2d (omega = k b U/L)",
           "source": [r["dir"] for r in cfg["runs"]]}
    json.dump(out, open(a.out, "w"), indent=1)
    print("%-10s %10s %10s %21s   (bootstrap over %s)" % ("constant", "fit", "sd", "10-90 %", a.boot_mode))
    for k in free:
        print("%-10s %10.4f %10s %21s" % (k, p[k], "%.4f" % sd[k] if sd[k] is not None else "-",
                                         "%.4f - %.4f" % tuple(ci[k]) if ci[k] is not None else "-"))
    print("power: c3 %.4f, c1 %.4f (rms %.4f W over %d beats)" % (pw["c3"], pw["c1"], prms, npw))
    for x in per:
        print("%-8s beats %2d  rms U %.4f m/s  turn %.4f rad/s  slip %.4f m/s" % (x["run"], x["beats"], x["rms_U"] or 0,
              x["rms_W"] or 0, x["rms_V"] or 0))
    if a.table:
        with open(a.table, "w") as fo:
            fo.write("run,t0,t1,f,a,U,W,V,P,use,U_m,W_m,V_m\n")
            for R in runs:
                for b, mm in zip(R["beats"], model_beats(R, p)):
                    fo.write("%s,%.4f,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f,%.5f,%d,%.5f,%.5f,%.5f\n" % (
                        R["name"], b["t0"], b["t1"], b["f"], b["a"], b["U"], b["W"], b["V"], b["P"], int(b["use"]), mm[0], mm[1], mm[2]))
    if a.ctrl:
        write_ctrl(p, a.ctrl)


def write_ctrl(p, path, body_id=0):
    """the settings of the brain of REEF3D for this body (X 200 s0 k, X 202 n_u n_d, X 203 sd, X 204 tau_c, X 209 C_d).
    With the drag law of the glide (C_d), X 209 can give it to the brain (patch 0041; written commented), and X 204
    is the exponential equivalent, L / (C_d U1 ln 2) with U1 the speed of the gait at f0: the same half-time of the
    glide from U1 (the brain without X 209 glides with it; with X 209 it is the decay of the water that passed the
    body at the start)"""
    s1 = bd.stride(p, p["f0"], 0.0)
    n_down, tau_c = p["n_down"], p["tau_coast"]
    # (comment lines start with "/", as REEF3D's ctrl.txt reads them)
    lines = ["// the brain's settings for the body of %s (fit_body.py)" % os.path.basename(path)]
    if p.get("C_d", 0.0) > 0.0:
        U1 = s1 * p["f0"] * p["L"]
        tau_c = p["L"] / (p["C_d"] * U1 * math.log(2.0))
        lines.append("// the glide by the drag law, C_d %.4f (X 209, patch 0041); X 204 its exponential equivalent from "
                     "U1 = %.3f m/s" % (p["C_d"], U1))
    lines += ["X 200 %d %.4f %.4f 0" % (body_id, s1, p["k"]),
              "X 202 %d %.3f %.3f" % (body_id, p["n_up"], n_down),
              "X 203 %d 0.05" % body_id,
              "X 204 %d %.2f" % (body_id, tau_c)]
    if p.get("C_d", 0.0) > 0.0:
        # the drag law of the glide for the brain: off (commented) by default; with the 2D swimmer it predicts the
        # glide better but holds the station no better (validation 31)
        lines.append("// X 209 %d %.4f" % (body_id, p["C_d"]))
    open(path, "w").write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
