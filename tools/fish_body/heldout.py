#!/usr/bin/env python3
# Plan 2, work package E: how well the agent fish stands in for the resolved fish. Each held-out resolved run with
# the brain (not in the fit) is repeated by an agent fish with the same brain and the fitted body (fit_body.py), and
# with each body of the bootstrap of the fit for the spread; the track, the distance to the target, the frequency of
# the gait and the time to the station are compared.
#
#   python3 heldout.py heldout.json body.json [--out heldout.json] [--csv tracks.csv] [--nboot 20]
#
# heldout.json: {"runs": [{"name": "cur", "dir": "<folder with the REEF3D_6DOF_* logs>", "T": 20,
#                          "flow": [0.2, 0.0], "walls": [0, 0, 1, 1], "target": [3.5, 3.0] | "targets": [[t, x, z], ...],
#                          "start": [1.6, 3.0, 0.0], "brain": {"U_des": 0.5, "b_max": 1.0, "D": 1.5}}, ...],
#                "box": [0, 8, 0, 6], "probes": 1.0, "compare": {"old": {<the body of fish_agents.py>}}}
#
# Architect: Hans Bihs
import argparse, json, os, sys, copy
import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..", "fish_agents"))
import fish_agents as fa                                             # noqa: E402


def table(p):
    rows = []
    for l in open(p):
        try:
            rows.append([float(x) for x in l.split()])
        except ValueError:
            pass
    n = max(set(len(r) for r in rows), key=lambda k: sum(len(r) == k for r in rows))
    return np.array([r for r in rows if len(r) == n])


def target_at(run, t):
    if "targets" in run:
        T = np.array(run["targets"], float)
        return np.interp(t, T[:, 0], T[:, 1]), np.interp(t, T[:, 0], T[:, 2])
    return np.full_like(t, run["target"][0]), np.full_like(t, run["target"][1])


def agent(run, body, box, probes):
    x, z, h = run.get("start", [1.6, 3.0, 0.0])
    case = {"T": run["T"], "dt": 0.01, "out_dt": 0.1, "box": box, "walls": run.get("walls", [1, 1, 1, 1]),
            "flow": {"type": "uniform", "u": run.get("flow", [0.0, 0.0])}, "brain": run.get("brain", {}),
            "fish": [{"x": x, "z": z, "h": h, "target": run.get("target", [x, z])}],
            "body": dict(body, probes=probes)}
    if "targets" in run:
        case["targets"] = run["targets"]
    _, rows = fa.run(case)
    R = np.array([r[:9] for r in rows], float)
    return R          # t id x z heading U f a b


def metrics(run, R, P, Bl):
    t = R[:, 0]
    xc, zc = np.interp(t, P[:, 0], P[:, 1]), np.interp(t, P[:, 0], P[:, 3])
    tx, tz = target_at(run, t)
    da, dc = np.hypot(R[:, 2] - tx, R[:, 3] - tz), np.hypot(xc - tx, zc - tz)
    err = np.hypot(R[:, 2] - xc, R[:, 3] - zc)
    m = t > 3.0
    # the time to the station: first within 0.25 L (station cases)
    def first(d):
        k = np.where(d < 0.25)[0]
        return float(t[k[0]]) if len(k) else float("nan")
    fc = np.interp(t, Bl[:, 0], Bl[:, 7])
    return {"track_mean": float(err[m].mean()), "track_max": float(err[m].max()), "track_end": float(err[-1]),
            "dist_mean_agent": float(da[m].mean()), "dist_mean_cfd": float(dc[m].mean()),
            "t_station_agent": first(da), "t_station_cfd": first(dc),
            "f_mean_agent": float(R[m, 6].mean()), "f_mean_cfd": float(fc[m].mean()),
            "f_rms_diff": float(np.sqrt(np.mean((R[m, 6] - fc[m]) ** 2)))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("heldout"); ap.add_argument("body"); ap.add_argument("--out", default="heldout_result.json")
    ap.add_argument("--csv", default=""); ap.add_argument("--nboot", type=int, default=20)
    a = ap.parse_args()
    cfg = json.load(open(a.heldout)); B = json.load(open(a.body))
    box, probes = cfg.get("box", [0, 8, 0, 6]), cfg.get("probes", 1.0)
    fitted = dict(B["body"]); fitted["model"] = 2
    boots = []
    for q in B.get("boot", [])[:a.nboot]:
        bb = dict(fitted); bb.update(zip(B["free"], q)); boots.append(bb)
    out, tracks = {}, []
    for run in cfg["runs"]:
        P = np.loadtxt(os.path.join(run["dir"], "REEF3D_6DOF_position_0.dat"), skiprows=1)
        _, ip = np.unique(P[:, 0], return_index=True); P = P[ip]
        Bl = table(os.path.join(run["dir"], "REEF3D_6DOF_brain_0.dat")); Bl = Bl[Bl[:, 0] > 0]
        res = {}
        R = agent(run, fitted, box, probes); res["fitted"] = metrics(run, R, P, Bl)
        tracks += [(run["name"], "fitted", r) for r in R]
        for lab, body in cfg.get("compare", {}).items():
            Rc = agent(run, body, box, probes if body.get("model", 1) == 2 else 0.0); res[lab] = metrics(run, Rc, P, Bl)
            tracks += [(run["name"], lab, r) for r in Rc]
        if boots:
            ms = [metrics(run, agent(run, bb, box, probes), P, Bl) for bb in boots]
            res["boot"] = {k: [float(np.nanpercentile([m[k] for m in ms], 10)), float(np.nanpercentile([m[k] for m in ms], 90))]
                           for k in ms[0]}
        out[run["name"]] = res
        f = res["fitted"]
        print("%-7s track error after 3 s mean %.2f m (max %.2f); distance to the target agent %.2f / CFD %.2f m; "
              "station at %.1f / %.1f s; f agent %.2f / CFD %.2f Hz" % (run["name"], f["track_mean"], f["track_max"],
              f["dist_mean_agent"], f["dist_mean_cfd"], f["t_station_agent"], f["t_station_cfd"], f["f_mean_agent"], f["f_mean_cfd"]))
        for lab in cfg.get("compare", {}):
            g = res[lab]
            print("        %-6s track error mean %.2f m; distance %.2f m; station at %.1f s; f %.2f Hz" % (
                lab, g["track_mean"], g["dist_mean_agent"], g["t_station_agent"], g["f_mean_agent"]))
        if "boot" in res:
            bt = res["boot"]
            print("        bootstrap 10-90 %%: track error %.2f-%.2f m, distance %.2f-%.2f m, f %.2f-%.2f Hz" % (
                *bt["track_mean"], *bt["dist_mean_agent"], *bt["f_mean_agent"]))
    json.dump(out, open(a.out, "w"), indent=1)
    if a.csv:
        with open(a.csv, "w") as fo:
            fo.write("run,body,t,x,z,heading,U,f,a,b\n")
            for name, lab, r in tracks:
                fo.write("%s,%s,%.2f,%.4f,%.4f,%.2f,%.4f,%.4f,%.4f,%.4f\n" % (name, lab, r[0], r[2], r[3], r[4], r[5], r[6], r[7], r[8]))


if __name__ == "__main__":
    main()
