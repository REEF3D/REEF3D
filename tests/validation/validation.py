#!/usr/bin/env python3
"""
REEF3D validation suite: cases with an analytical (or reference) solution.

Unlike the regression suite (does the code give the same numbers?), these cases check that the
numbers are right. They use the same case format and runner as tests/regression (case.json with
control.txt/ctrl.txt or "base" + overrides), plus

  "time":  simulated time (instead of a number of steps)
  "check": {"type": <checker>, ...parameters}   see CHECKS below
  "check_set": {...}   override single check parameters of the base case

Usage
-----
  ./validation.py list
  ./validation.py run   --reef3d BIN --divemesh DM --out DIR [--cases 'channel_*'] [--tags ...]
  ./validation.py check DIR [--cases ...]          # evaluate an existing run
  ./validation.py run ... --check                  # run + evaluate

Result per case: PASS / FAIL with the error measure of the checker and its tolerance.
Pure Python 3 standard library.
"""

import argparse
import glob
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "regression"))
import regression as R  # noqa: E402

R.CASES_DIR = os.path.join(HERE, "cases")


# ----------------------------------------------------------------------------------------------
# probe output
# ----------------------------------------------------------------------------------------------

def read_probes(rundir):
    """REEF3D_CFD_ProbePoint/*-<n>.dat -> {n: {"xyz": (x,y,z), "t": [...], "u": [...], "v":..., "w":..., "p":...}}"""
    out = {}
    for fn in glob.glob(os.path.join(rundir, "REEF3D_CFD_ProbePoint", "*.dat")):
        n = int(os.path.splitext(fn)[0].rsplit("-", 1)[1])
        xyz = None
        rows = []
        with open(fn) as f:
            lines = f.read().split("\n")
        for i, l in enumerate(lines):
            t = l.split()
            if l.startswith("x_coord") and i + 1 < len(lines):
                xyz = tuple(float(x) for x in lines[i + 1].split()[1:4])
            if len(t) >= 5:
                try:
                    rows.append([float(x) for x in t])
                except ValueError:
                    pass
        out[n] = {"xyz": xyz, "t": [r[0] for r in rows], "u": [r[1] for r in rows],
                  "v": [r[2] for r in rows], "w": [r[3] for r in rows], "p": [r[4] for r in rows]}
    return out


def interp(ts, ys, t):
    if not ts:
        return float("nan")
    if t <= ts[0]:
        return ys[0]
    for i in range(1, len(ts)):
        if ts[i] >= t:
            a = (t - ts[i - 1]) / (ts[i] - ts[i - 1]) if ts[i] > ts[i - 1] else 0.0
            return ys[i - 1] + a * (ys[i] - ys[i - 1])
    return ys[-1] if abs(ts[-1] - t) < 1e-9 + 1e-6 * t else float("nan")


# ----------------------------------------------------------------------------------------------
# checkers
# ----------------------------------------------------------------------------------------------

def channel_startup_exact(z, t, g, nu, H, nmax=401):
    """Body-force driven laminar channel flow between no-slip walls at z=0 and z=H, from rest:
       u_t = g + nu u_zz.  u = g/(2nu) z(H-z) - sum_{n odd} 4 g H^2/(nu n^3 pi^3) sin(n pi z/H) exp(-n^2 pi^2 nu t/H^2)"""
    us = g / (2.0 * nu) * z * (H - z)
    s = 0.0
    for n in range(1, nmax + 1, 2):
        s += 4.0 * g * H * H / (nu * n ** 3 * math.pi ** 3) * math.sin(n * math.pi * z / H) \
            * math.exp(-n * n * math.pi * math.pi * nu * t / (H * H))
    return us - s


def check_channel_startup(c, rundir):
    """Probe velocities u(z,t) against channel_startup_exact at the given times.
    error = max over probes of |u - u_exact| / u_max_steady, separately for the last time
    ("final", close to steady state) and the earlier times ("transient").
    Tolerances: "tol_final", "tol_transient" (or one "tol" for both)."""
    ck = c["check"]
    g, nu, H, z0 = ck["g"], ck["nu"], ck.get("H", 1.0), ck.get("z0", 0.0)
    umax = g * H * H / (8.0 * nu)
    probes = read_probes(rundir)
    if not probes:
        return {"ok": False, "error": float("nan"), "note": "no probe output"}
    rows = []
    err_t = 0.0
    err_s = 0.0
    for n in sorted(probes):
        pr = probes[n]
        z = pr["xyz"][2] - z0
        for t in ck["times"]:
            u = interp(pr["t"], pr["u"], t)
            ue = channel_startup_exact(z, t, g, nu, H)
            e = abs(u - ue) / umax
            rows.append((n, z, t, u, ue, e))
            if t == ck["times"][-1]:
                err_s = max(err_s, e)
            else:
                err_t = max(err_t, e)
    err = max(err_t, err_s)
    tf = ck.get("tol_final", ck.get("tol"))
    tt = ck.get("tol_transient", ck.get("tol"))
    return {"ok": err_s <= tf and err_t <= tt, "error": err, "error_transient": err_t, "error_final": err_s,
            "tol": "final %g, transient %g" % (tf, tt), "rows": rows}


def check_conserved_integral(c, rundir):
    """Volume integral of a field (sum over cells of field*vol, from the regression state dumps)
    must stay constant: error = |I_end - I_start| / |I_start|. Needs the "vol" field in the dump."""
    ck = c["check"]
    sf = R.state_files(rundir)
    if len(sf) < 2:
        return {"ok": False, "error": float("nan"), "note": "need initial and final state dump"}
    def integral(count):
        tot = []
        for r in sorted(sf[count]):
            st = R.read_state(sf[count][r])
            f, v = st["fields"][ck["field"]], st["fields"]["vol"]
            tot.extend(x * y for x, y in zip(f, v))
        return math.fsum(tot)
    i0, i1 = integral(min(sf)), integral(max(sf))
    err = abs(i1 - i0) / max(abs(i0), 1e-300)
    return {"ok": err <= ck["tol"], "error": err, "tol": ck["tol"],
            "rows_text": "integral of %s: start %.12g, end %.12g" % (ck["field"], i0, i1)}


def check_max_abs(c, rundir):
    """Maximum absolute value of a dumped field in the final state (e.g. "div": the discrete
    divergence after the projection) must be below "tol"."""
    ck = c["check"]
    st = R.merged_final_state(rundir)
    if st is None:
        return {"ok": False, "error": float("nan"), "note": "no state dump"}
    m = 0.0
    for s in st[1]:
        f = s["fields"].get(ck["field"])
        if f is None:
            return {"ok": False, "error": float("nan"), "note": "field %s not in dump" % ck["field"]}
        m = max(m, max((abs(x) for x in f), default=0.0))
    return {"ok": m <= ck["tol"], "error": m, "tol": ck["tol"],
            "rows_text": "max |%s| = %.3e (final step %d)" % (ck["field"], m, st[0])}


CHECKS = {"channel_startup": check_channel_startup,
          "max_abs": check_max_abs,
          "conserved_integral": check_conserved_integral}


# ----------------------------------------------------------------------------------------------

def evaluate(names, outroot):
    results = []
    for n in names:
        c = R.load_case(n)
        rd = os.path.join(outroot, n)
        st = R.run_status(rd)
        if st.get("status") != "ok":
            results.append((n, {"ok": False, "error": float("nan"), "note": st.get("status")}))
            continue
        if c.get("check_set"):
            c["check"] = dict(c["check"], **c["check_set"])
        results.append((n, CHECKS[c["check"]["type"]](c, rd)))
    lines = ["# REEF3D validation", "", "run: `%s`" % outroot, "",
             "| case | result | error | transient | final | tol |", "|---|---|---|---|---|---|"]
    for n, r in results:
        lines.append("| %s | %s | %.3g | %s | %s | %s |" % (
            n, "PASS" if r["ok"] else "**FAIL**", r["error"],
            "%.3g" % r["error_transient"] if "error_transient" in r else "",
            "%.3g" % r["error_final"] if "error_final" in r else "", r.get("tol", r.get("note", ""))))
    lines += ["", "## Details", ""]
    for n, r in results:
        if "rows_text" in r:
            lines += ["### %s" % n, "", r["rows_text"], ""]
        if "rows" not in r:
            continue
        lines += ["### %s" % n, "", "| probe | z | t | u | u exact | err/umax |", "|---|---|---|---|---|---|"]
        for row in r["rows"]:
            lines.append("| %d | %.4g | %.4g | %.5g | %.5g | %.3g |" % row)
        lines.append("")
    txt = "\n".join(lines)
    with open(os.path.join(outroot, "validation.md"), "w") as f:
        f.write(txt + "\n")
    with open(os.path.join(outroot, "validation.json"), "w") as f:
        json.dump([{"case": n, **{k: v for k, v in r.items() if k != "rows"}} for n, r in results], f, indent=1)
    print(txt.split("\n## Details")[0])
    print("\nreport: %s" % os.path.join(outroot, "validation.md"))
    return 0 if all(r["ok"] for _, r in results) else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("list"); p.add_argument("--cases", nargs="*"); p.add_argument("--tags", nargs="*")
    p = sub.add_parser("run"); p.add_argument("--cases", nargs="*"); p.add_argument("--tags", nargs="*")
    p.add_argument("--reef3d", required=True); p.add_argument("--divemesh", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--mpirun", default=os.environ.get("REEF3D_MPIRUN", "mpirun"))
    p.add_argument("--timeout", type=int, default=3600); p.add_argument("--keep", action="store_true")
    p.add_argument("--check", action="store_true", help="evaluate after running")
    p = sub.add_parser("check"); p.add_argument("run"); p.add_argument("--cases", nargs="*")
    p.add_argument("--tags", nargs="*")
    a = ap.parse_args()

    if a.cmd == "list":
        a.steps = None
        return R.cmd_list(a)
    if a.cmd == "run":
        a.steps = None
        rc = R.cmd_run(a)
        if a.check:
            rc |= evaluate(R.select_cases(a.cases, a.tags), a.out)
        return rc
    if a.cmd == "check":
        names = [n for n in R.select_cases(a.cases, a.tags) if os.path.isdir(os.path.join(a.run, n))]
        return evaluate(names, a.run)


if __name__ == "__main__":
    sys.exit(main())
