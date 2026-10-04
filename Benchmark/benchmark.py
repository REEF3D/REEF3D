#!/usr/bin/env python3
"""
REEF3D literature benchmark suite.

Established cases from the literature (laboratory experiments, analytical and semi-analytical
solutions) for REEF3D::CFD, ::NHFLOW, ::FNPF and ::SFLOW. Each case is run on two levels:

  nightly   coarse grid / short run, minutes on 1-4 ranks, loose tolerances
  release   the published (paper/tutorial) resolution, hours on a cluster, tight tolerances

The case format follows Regression/ (case.json + control.txt + ctrl.txt,
or "base" + key overrides), plus

  "levels":    {"nightly": {...}, "release": {...}}   per-level overrides of any case key
               (control_set, ctrl_set, control_add, ctrl_add, np, time, check_set, files)
  "check":     {"type": <checker>, ...}              how to evaluate the run (CHECKS below)
  "reference": {"citation": ..., "data": [...], "notes": ...}   literature source

Usage
-----
  ./benchmark.py list   [--level nightly] [--cases 'cfd_*'] [--tags ...]
  ./benchmark.py write  --out DIR [--level nightly] [--cases ...]        resolved inputs only
  ./benchmark.py run    --reef3d BIN --divemesh DM --out DIR [--level nightly] [--cases ...] [--check]
  ./benchmark.py check  DIR [--cases ...]                                evaluate an existing run
  ./benchmark.py plot   DIR [--cases ...]                                PNG per case (needs matplotlib)

The level is stored in DIR/suite.json, so check/plot use the tolerances of the level that was run.
Pure Python 3 standard library; matplotlib only for plots.
"""

import argparse
import bisect
import fnmatch
import glob
import json
import math
import os
import shlex
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
CASES_DIR = os.path.join(HERE, "cases")
REFDATA_DIR = os.path.join(HERE, "refdata")
LEVELS = ("nightly", "release")
MODEL = {2: "SFLOW", 3: "FNPF", 5: "NHFLOW", 6: "CFD"}

# bulk output removed after a run unless --keep
BULK_OUTPUT = ["REEF3D_CFD_VTU", "REEF3D_CFD_6DOF_VTP", "REEF3D_CFD_6DOF_Normals_VTP", "REEF3D_CFD_FSF",
               "REEF3D_CFD_Topo", "REEF3D_CFD_Sediment_VTP",
               "REEF3D_NHFLOW_VTU", "REEF3D_NHFLOW_VTP_FSF", "REEF3D_NHFLOW_VTP_BED",
               "REEF3D_FNPF_VTU", "REEF3D_FNPF_VTP_FSF", "REEF3D_FNPF_VTP_BED",
               "REEF3D_SFLOW_VTP", "REEF3D_SFLOW_VTP_FSF", "REEF3D_SFLOW_VTP_BED", "DIVEMesh_Paraview"]


# ==============================================================================================
# cases
# ==============================================================================================

def key_of(line):
    t = line.split()
    if len(t) >= 2 and len(t[0]) == 1 and t[0].isalpha() and t[1].isdigit():
        return t[0] + " " + t[1]
    return None


def apply_overrides(lines, setmap, addlist):
    """setmap {"N 40": "13"} replaces every line with that key (null removes, a list writes
    several lines); addlist lines are appended as they are."""
    out = [l for l in lines if key_of(l) not in setmap]
    for k, v in setmap.items():
        if v is None:
            continue
        for x in (v if isinstance(v, list) else [v]):
            out.append("%s %s" % (k, x))
    return out + list(addlist)


def read_lines(path):
    with open(path) as f:
        return [l.rstrip("\n") for l in f]


def merge_level(c, level):
    """apply the overrides of one level on top of the case keys"""
    lv = (c.get("levels") or {}).get(level)
    if lv is None:
        return None
    c = dict(c)
    for k, v in lv.items():
        if k in ("control_set", "ctrl_set", "check_set"):
            c[k] = dict(c.get(k, {}), **v)
        elif k in ("control_add", "ctrl_add", "files"):
            c[k] = list(c.get(k, [])) + list(v)
        else:
            c[k] = v
    return c


def load_case(name, level="release", _seen=None):
    """case with 'base' inheritance and the level overrides resolved; None if the case has no such level"""
    _seen = _seen or set()
    if name in _seen:
        raise ValueError("circular base in case %s" % name)
    _seen.add(name)
    cdir = os.path.join(CASES_DIR, name)
    with open(os.path.join(cdir, "case.json")) as f:
        raw = json.load(f)
    if raw.get("base"):
        b = load_case(raw["base"], None, _seen)
        control, ctrl, files = list(b["control_lines"]), list(b["ctrl_lines"]), list(b["files_abs"])
        c = dict(raw)
        for k, v in b.items():
            if k not in ("name", "dir", "base", "control_lines", "ctrl_lines", "files_abs", "files",
                         "control_set", "ctrl_set", "control_add", "ctrl_add", "description", "check_set"):
                c.setdefault(k, v)
        if "check" in raw and "check" in b and raw["check"].get("type") in (None, b["check"]["type"]):
            c["check"] = dict(b["check"], **raw["check"])
    else:
        c = dict(raw)
        control = read_lines(os.path.join(cdir, "control.txt"))
        ctrl = read_lines(os.path.join(cdir, "ctrl.txt"))
        files = []
    c["name"], c["dir"] = name, cdir
    if level is not None:
        c = merge_level(c, level)
        if c is None:
            return None
        c["level"] = level
    c["control_lines"] = apply_overrides(control, c.get("control_set", {}), c.get("control_add", []))
    c["ctrl_lines"] = apply_overrides(ctrl, c.get("ctrl_set", {}), c.get("ctrl_add", []))
    c["files_abs"] = files + [os.path.join(cdir, x) for x in c.get("files", [])]
    if c.get("check_set"):
        c["check"] = dict(c.get("check", {}), **c["check_set"])
    c.setdefault("np", 1)
    c.setdefault("tags", [])
    return c


def all_case_names():
    return sorted(d for d in os.listdir(CASES_DIR) if os.path.isfile(os.path.join(CASES_DIR, d, "case.json")))


def select_cases(patterns, tags, level):
    out = []
    for n in all_case_names():
        if patterns and not any(fnmatch.fnmatch(n, p) for p in patterns):
            continue
        c = load_case(n, level)
        if c is None or c.get("abstract"):
            continue
        if tags and not set(tags) & set(c["tags"]):
            continue
        out.append(n)
    return out


def ctrl_value(c, key, idx=0, default=None):
    for l in c["ctrl_lines"]:
        if key_of(l) == key:
            t = l.split()
            return float(t[2 + idx])
    return default


def model_of(c):
    return MODEL[int(ctrl_value(c, "A 10", 0, 6))]


# ==============================================================================================
# running
# ==============================================================================================

PRINT_KEYS_REMOVED = ["P 10", "P 20", "P 30", "P 35", "P 40", "P 41", "P 42", "P 180", "P 182", "P 185"]


def write_inputs(c, rundir, keep_vtu=False):
    os.makedirs(rundir, exist_ok=True)
    np_ = str(c["np"])
    control = apply_overrides(c["control_lines"], {"M 10": np_}, [])
    sets = {"M 10": np_}
    if c.get("time") is not None:
        sets["N 41"] = repr(float(c["time"]))
    if not keep_vtu:
        for k in PRINT_KEYS_REMOVED:
            sets[k] = None
    ctrl = apply_overrides(c["ctrl_lines"], sets, [])
    with open(os.path.join(rundir, "control.txt"), "w") as f:
        f.write("\n".join(control) + "\n")
    with open(os.path.join(rundir, "ctrl.txt"), "w") as f:
        f.write("\n".join(ctrl) + "\n")
    for src in c["files_abs"]:
        shutil.copy(src, rundir)
    gen = c.get("generate")
    if gen:  # input files made by a tool script (e.g. geo.dat bathymetry)
        cmd = [sys.executable, os.path.join(HERE, "tools", gen["script"])] + [str(a) for a in gen.get("args", [])]
        subprocess.run(cmd, cwd=rundir, check=True)


def run_case(c, args, outroot):
    rundir = os.path.join(outroot, c["name"])
    if os.path.isdir(rundir):
        shutil.rmtree(rundir)
    write_inputs(c, rundir, args.keep)
    info = {"case": c["name"], "level": c["level"], "np": c["np"], "status": "ok",
            "reef3d": os.path.abspath(args.reef3d), "divemesh": os.path.abspath(args.divemesh),
            "time_divemesh": 0.0, "time_reef3d": 0.0}
    env = dict(os.environ)
    env.setdefault("OMP_NUM_THREADS", "1")
    t0 = time.time()
    with open(os.path.join(rundir, "divemesh.log"), "w") as log:
        r = subprocess.run([os.path.abspath(args.divemesh)], cwd=rundir, stdout=log, stderr=subprocess.STDOUT)
    info["time_divemesh"] = round(time.time() - t0, 2)
    if r.returncode != 0:
        info["status"] = "divemesh failed (%d)" % r.returncode
        return finish_run(rundir, info, args)
    cmd = shlex.split(args.mpirun) + ["-np", str(c["np"]), os.path.abspath(args.reef3d)]
    t0 = time.time()
    try:
        with open(os.path.join(rundir, "reef3d.log"), "w") as log:
            r = subprocess.run(cmd, cwd=rundir, stdout=log, stderr=subprocess.STDOUT,
                               timeout=args.timeout, env=env)
        if r.returncode != 0:
            info["status"] = "reef3d failed (%d)" % r.returncode
    except subprocess.TimeoutExpired:
        info["status"] = "reef3d timeout (%ds)" % args.timeout
    info["time_reef3d"] = round(time.time() - t0, 2)
    return finish_run(rundir, info, args)


def finish_run(rundir, info, args):
    if not args.keep:
        for d in BULK_OUTPUT:
            shutil.rmtree(os.path.join(rundir, d), ignore_errors=True)
        for d in glob.glob(os.path.join(rundir, "DIVEMesh_Grid*")):
            shutil.rmtree(d, ignore_errors=True)
        for f in glob.glob(os.path.join(rundir, "grid-*.dat")):
            os.remove(f)
    with open(os.path.join(rundir, "run.json"), "w") as f:
        json.dump(info, f, indent=1)
    return info


def run_status(d):
    try:
        with open(os.path.join(d, "run.json")) as f:
            return json.load(f)
    except OSError:
        return {"status": "not run"}


# ==============================================================================================
# reading REEF3D output and reference data
# ==============================================================================================

def is_num(x):
    try:
        float(x)
        return True
    except ValueError:
        return False


def numeric_rows(path, ncol=None):
    rows = []
    with open(path, errors="replace") as f:
        for l in f:
            t = l.split()
            if t and all(is_num(x) for x in t) and (ncol is None or len(t) == ncol):
                rows.append([float(x) for x in t])
    return rows


def read_refdata(rel):
    """reference data file in refdata/: numeric rows (comment lines start with #)"""
    return numeric_rows(os.path.join(REFDATA_DIR, rel))


def read_gauges(rundir, model, theory=False):
    """REEF3D_<M>_WSF/REEF3D-<M>-WSF-HG[-THEORY].dat -> (xs, ys, t, [eta per gauge]) or None"""
    fn = os.path.join(rundir, "REEF3D_%s_WSF" % model,
                      "REEF3D-%s-WSF-HG%s.dat" % (model, "-THEORY" if theory else ""))
    if not os.path.exists(fn):
        return None
    with open(fn) as f:
        lines = f.read().split("\n")
    ng = int(lines[0].split(":")[1])
    xs, ys, ts, cols = [], [], [], [[] for _ in range(ng)]
    i = 1
    while i < len(lines) and len(xs) < ng:
        t = lines[i].split()
        if len(t) == 3 and all(is_num(x) for x in t):
            xs.append(float(t[1]))
            ys.append(float(t[2]))
        i += 1
    for l in lines[i:]:
        t = l.split()
        if len(t) == ng + 1 and all(is_num(x) for x in t):
            ts.append(float(t[0]))
            for g in range(ng):
                cols[g].append(float(t[g + 1]))
    return xs, ys, ts, cols


def read_probe_file(path):
    """probe files: header lines, then 'time value(s)' rows -> (t, [col1, col2, ...])"""
    rows = [r for r in numeric_rows(path) if len(r) >= 2]
    # drop the coordinate line (id x y z) written before the data
    rows = [r for r in rows if len(r) == len(rows[-1])] if rows else rows
    if rows and len(rows) > 1 and rows[0][0] > rows[1][0]:
        rows = rows[1:]
    t = [r[0] for r in rows]
    cols = [[r[k] for r in rows] for k in range(1, len(rows[0]))] if rows else []
    return t, cols


def read_force(rundir, model, n=1):
    """REEF3D_<M>_Force/REEF3D_<M>_Force-n.dat: it, time, Fx, Fy, Fz"""
    fn = os.path.join(rundir, "REEF3D_%s_Force" % model, "REEF3D_%s_Force-%d.dat" % (model, n))
    if not os.path.exists(fn):
        return None
    rows = [r for r in numeric_rows(fn) if len(r) == 5]
    return {"t": [r[1] for r in rows], "Fx": [r[2] for r in rows], "Fy": [r[3] for r in rows],
            "Fz": [r[4] for r in rows]}


def read_6dof_position(rundir, model, n=0):
    fns = sorted(glob.glob(os.path.join(rundir, "REEF3D_%s_6DOF" % model, "REEF3D_6DOF_position_%d.dat" % n)))
    if not fns:
        return None
    rows = [r for r in numeric_rows(fns[0]) if len(r) == 7]
    return {k: [r[i] for r in rows] for i, k in enumerate(["t", "x", "y", "z", "phi", "theta", "psi"])}


def read_wsflines(rundir, model):
    """REEF3D_<M>_WSFLINE/*.dat -> sorted list of (simtime, x[], eta[]) for the first line"""
    out = []
    for fn in glob.glob(os.path.join(rundir, "REEF3D_%s_WSFLINE" % model, "*.dat")):
        simtime = None
        xs, es = [], []
        with open(fn) as f:
            for l in f:
                if l.startswith("simtime"):
                    simtime = float(l.split(":")[1])
                    continue
                t = l.split()
                if len(t) >= 2 and is_num(t[0]) and is_num(t[1]) and simtime is not None:
                    if l.count("\t") >= 1 and len(t) <= 3:
                        xs.append(float(t[0]))
                        es.append(float(t[1]))
        if simtime is not None and xs:
            # the header block "1  y" also matches; keep rows from the data block only (x sorted)
            pairs = sorted(zip(xs, es))
            out.append((simtime, [p[0] for p in pairs], [p[1] for p in pairs]))
    out.sort(key=lambda r: r[0])
    return out


# ==============================================================================================
# numerics helpers
# ==============================================================================================

def interp(ts, ys, t):
    """linear interpolation in a sorted series; nan outside"""
    if not ts or t < ts[0] or t > ts[-1]:
        return float("nan")
    i = bisect.bisect_left(ts, t)
    if i == 0:
        return ys[0]
    a = (t - ts[i - 1]) / (ts[i] - ts[i - 1]) if ts[i] > ts[i - 1] else 0.0
    return ys[i - 1] + a * (ys[i] - ys[i - 1])


def rms(v):
    v = [x for x in v if x == x]
    return math.sqrt(sum(x * x for x in v) / len(v)) if v else float("nan")


def window(ts, ys, t0, t1):
    return [y for t, y in zip(ts, ys) if t0 <= t <= t1]


def zero_up_crossings(ts, ys, mean=None):
    if mean is None:
        mean = sum(ys) / len(ys)
    out = []
    for i in range(1, len(ys)):
        a, b = ys[i - 1] - mean, ys[i] - mean
        if a < 0.0 <= b and b != a:
            out.append(ts[i - 1] + (ts[i] - ts[i - 1]) * (-a) / (b - a))
    return out


def wave_height(ts, ys, t0, t1):
    """mean zero-up-crossing wave height in the window (max - min if fewer than 2 crossings)"""
    tt = [t for t in ts if t0 <= t <= t1]
    yy = window(ts, ys, t0, t1)
    if len(yy) < 3:
        return float("nan")
    zc = zero_up_crossings(tt, yy)
    if len(zc) < 2:
        return max(yy) - min(yy)
    hs = []
    for a, b in zip(zc[:-1], zc[1:]):
        seg = [y for t, y in zip(tt, yy) if a <= t <= b]
        if seg:
            hs.append(max(seg) - min(seg))
    return sum(hs) / len(hs)


def extrema(ts, ys):
    """local maxima and minima (t, y), parabolic refinement"""
    out = []
    for i in range(1, len(ys) - 1):
        if (ys[i] > ys[i - 1] and ys[i] >= ys[i + 1]) or (ys[i] < ys[i - 1] and ys[i] <= ys[i + 1]):
            y0, y1, y2 = ys[i - 1], ys[i], ys[i + 1]
            d = y0 - 2 * y1 + y2
            s = 0.5 * (y0 - y2) / d if d != 0 else 0.0
            s = max(-1.0, min(1.0, s))
            h = 0.5 * (ts[i + 1] - ts[i - 1])
            out.append((ts[i] + s * h, y1 - 0.25 * (y0 - y2) * s, "max" if y1 > y0 else "min"))
    return out


def gauge_index(g, spec):
    """spec {"index": n} (1-based, order of P 51 in ctrl.txt) or {"x": .., "y": ..} (nearest gauge)"""
    xs, ys = g[0], g[1]
    if "index" in spec:
        return spec["index"] - 1
    best, bi = 1e30, None
    for i, (x, y) in enumerate(zip(xs, ys)):
        d = (x - spec["x"]) ** 2 + (y - spec.get("y", y)) ** 2
        if d < best:
            best, bi = d, i
    return bi


def eta_series(c, g, i):
    """free-surface elevation at gauge i relative to the still water level: the value at the first
    output time is subtracted (all wave cases start from still water)"""
    ys = g[3][i]
    datum = c["check"].get("datum", "initial")
    d = ys[0] if datum == "initial" else float(datum)
    return [y - d for y in ys]


def result(ok, error, tol, text, **kw):
    r = {"ok": bool(ok), "error": error, "tol": tol, "rows_text": text}
    r.update(kw)
    return r


def fail(note):
    return {"ok": False, "error": float("nan"), "tol": "", "note": note, "rows_text": note}


# ==============================================================================================
# checkers
# ==============================================================================================

def sim_signal(c, rundir, s):
    """simulated signal for a 'signals' entry: (t[], y[]) in SI units"""
    model = model_of(c)
    src = s.get("source", "gauge")
    if src == "gauge":
        g = read_gauges(rundir, model)
        if g is None:
            return None
        i = gauge_index(g, s)
        return g[2], eta_series(c, g, i)
    if src == "force":
        f = read_force(rundir, model, s.get("n", 1))
        if f is None:
            return None
        k = s.get("component", "Fx")
        return f["t"], [s.get("factor", 1.0) * v for v in f[k]]
    if src == "pressure":
        if model == "CFD":
            fn = os.path.join(rundir, "REEF3D_CFD_PressureProbe", "REEF3D-CFD-Probe-Pressure-%d.dat" % s["n"])
        else:
            fn = os.path.join(rundir, "REEF3D_NHFLOW_ProbePoint", "REEF3D-NHFLOW-Probe_Press-%d.dat" % s["n"])
        if not os.path.exists(fn):
            return None
        t, cols = read_probe_file(fn)
        p0 = cols[0][0] if s.get("subtract_initial") else 0.0
        return t, [s.get("factor", 1.0) * (v - p0) for v in cols[0]]
    raise ValueError("unknown signal source %s" % src)


def ref_signal(s, window=None):
    """measured series in SI units (data_t_factor/offset, data_y_factor), optionally limited to a
    window [t0, t1] (in SI time of the data)"""
    rows = read_refdata(s["data"])
    tf, to, yf = s.get("data_t_factor", 1.0), s.get("data_t_offset", 0.0), s.get("data_y_factor", 1.0)
    rows.sort(key=lambda r: r[0])
    t, y = [r[0] * tf + to for r in rows], [r[1] * yf for r in rows]
    if window:
        keep = [i for i, x in enumerate(t) if window[0] <= x <= window[1]]
        t, y = [t[i] for i in keep], [y[i] for i in keep]
    return t, y


def check_timeseries(c, rundir):
    """Simulated time series (wave gauges, forces, pressure probes) against measured series.
    The measured time origin is arbitrary for the wave-tank data sets, so a common time lag is
    fitted (least squares on the 'align' signals, lag in 'lag_range'); for transient cases with a
    physical time origin (dam break) use "align": [] and lag 0.
    Optionally a small extra lag per signal ("local_lag", +-s) absorbs the uncertainty of gauge
    positions and of the digitised time axis of the measured series.
    Per signal: relative rms error  e = rms(sim - data) / rms(data)  over the measured window,
    and the height ratio  (max - min)_sim / (max - min)_data.  Pass if e <= tol_rms and
    |ratio - 1| <= tol_height for every signal."""
    ck = c["check"]
    sigs = []
    for s in ck["signals"]:
        sim = sim_signal(c, rundir, s)
        if sim is None or not sim[0]:
            return fail("no simulated output for %s" % s.get("name"))
        sigs.append((s, sim, ref_signal(s, ck.get("window"))))

    def err_for(lag, subset):
        tot, n = 0.0, 0
        for s, sim, ref in subset:
            for t, y in zip(*ref):
                v = interp(sim[0], sim[1], t + lag)
                if v == v:
                    tot += (v - y) ** 2
                    n += 1
        return tot / n if n else float("inf")

    align = [x for x in sigs if x[0].get("name") in ck.get("align", [])]
    lag = ck.get("lag", 0.0)
    if align:
        lo, hi = ck["lag_range"]
        step = ck.get("lag_step", 0.005)
        best = (float("inf"), 0.0)
        k = 0
        while lo + k * step <= hi:
            L = lo + k * step
            e = err_for(L, align)
            if e < best[0]:
                best = (e, L)
            k += 1
        lag = best[1]
    ok, worst, txt = True, 0.0, ["time lag sim - data: %.3f s" % lag, ""]
    txt.append("| signal | rel. rms error | tol | height ratio | tol | extra lag [s] | |")
    txt.append("|---|---|---|---|---|---|---|")
    per = {}
    for s, sim, ref in sigs:
        # optional small per-signal lag on top of the common one (gauge position / digitising
        # uncertainty of the measured series), searched in +-local_lag
        ll = s.get("local_lag", ck.get("local_lag", 0.0))
        dl = 0.0
        if ll > 0:
            cand = [k * 0.002 for k in range(-int(ll / 0.002), int(ll / 0.002) + 1)]
            dl = min(cand, key=lambda d: err_for(lag + d, [(s, sim, ref)]))
        diff, yy, ss = [], [], []
        for t, y in zip(*ref):
            v = interp(sim[0], sim[1], t + lag + dl)
            if v == v:
                diff.append(v - y)
                yy.append(y)
                ss.append(v)
        if len(diff) < max(3, 0.8 * len(ref[0])):
            return fail("simulation does not cover the measured window for %s" % s.get("name"))
        e = rms(diff) / rms([y - sum(yy) / len(yy) for y in yy]) if ck.get("rms_about_mean", True) else rms(diff) / rms(yy)
        hr = (max(ss) - min(ss)) / (max(yy) - min(yy))
        te = s.get("tol_rms", ck.get("tol_rms"))
        th = s.get("tol_height", ck.get("tol_height"))
        good = (te is None or e <= te) and (th is None or abs(hr - 1) <= th)
        ok &= good
        worst = max(worst, e)
        per[s.get("name")] = {"rms": e, "height_ratio": hr, "local_lag": dl}
        txt.append("| %s | %.3f | %s | %.3f | %s | %+.3f | %s |" % (s.get("name"), e, te, hr,
                                                                "±%s" % th if th is not None else "-", dl,
                                                                "" if good else "**FAIL**"))
    if "peak" in ck:  # transient: peak value and its time for selected signals
        for s, sim, ref in sigs:
            if s.get("name") not in ck["peak"]["signals"]:
                continue
            t0, t1 = ck["peak"]["window"]
            iy = max(range(len(ref[1])), key=lambda i: ref[1][i] if t0 <= ref[0][i] <= t1 else -1e30)
            iz = max(range(len(sim[1])), key=lambda i: sim[1][i] if t0 <= sim[0][i] - lag <= t1 else -1e30)
            pr = sim[1][iz] / ref[1][iy]
            dt = (sim[0][iz] - lag) - ref[0][iy]
            good = abs(pr - 1) <= ck["peak"]["tol_value"] and abs(dt) <= ck["peak"]["tol_time"]
            ok &= good
            txt.append("")
            txt.append("%s peak: sim %.4g at %.3f s, data %.4g at %.3f s -> ratio %.3f (tol ±%g), dt %.3f s (tol ±%g) %s" % (
                s.get("name"), sim[1][iz], sim[0][iz] - lag, ref[1][iy], ref[0][iy], pr, ck["peak"]["tol_value"], dt,
                ck["peak"]["tol_time"], "" if good else "**FAIL**"))
    return result(ok, worst, "rms %s / height %s" % (ck.get("tol_rms"), ck.get("tol_height")), "\n".join(txt),
                  lag=lag, signals=per)


def check_theory_gauges(c, rundir):
    """P 51 against P 50 (the wave theory used for generation) at the same positions: the wave
    should propagate along the tank without losing height or phase (e.g. Fenton's 5th-order Stokes
    theory). Error per gauge: rms(eta - eta_theory)/rms(eta_theory) and height ratio in the window."""
    ck = c["check"]
    model = model_of(c)
    g, th = read_gauges(rundir, model), read_gauges(rundir, model, theory=True)
    if g is None or th is None:
        return fail("no gauge or theory output (P 51 / P 50)")
    t0, t1 = ck["t_start"], ck["t_end"]
    ok, worst = True, 0.0
    txt = ["| x | rel. rms error | height ratio | |", "|---|---|---|---|"]
    for i, x in enumerate(g[0]):
        j = min(range(len(th[0])), key=lambda k: abs(th[0][k] - x))
        sim = eta_series(c, g, i)
        ref = th[3][j]
        tt = [t for t in g[2] if t0 <= t <= t1]
        d = [interp(g[2], sim, t) - interp(th[2], ref, t) for t in tt]
        r = [interp(th[2], ref, t) for t in tt]
        e = rms(d) / rms(r)
        hr = wave_height(g[2], sim, t0, t1) / wave_height(th[2], ref, t0, t1)
        good = e <= ck["tol_rms"] and abs(hr - 1) <= ck["tol_height"]
        ok &= good
        worst = max(worst, e)
        txt.append("| %.2f | %.3f | %.3f | %s |" % (x, e, hr, "" if good else "**FAIL**"))
    return result(ok, worst, "rms %g / height ±%g" % (ck["tol_rms"], ck["tol_height"]), "\n".join(txt))


def check_front_arrival(c, rundir):
    """Dam break surge front (Martin & Moyce 1952): wave gauges at x = x_wall + Z a; the front has
    arrived when the surface at the gauge is h_thr above the bed. The arrival time, as
    T = t sqrt(2g/a), is compared with the measured T(Z) (linear interpolation of the data).
    Error: max over gauges of |T_sim - T_exp| / T_exp (and the mean)."""
    ck = c["check"]
    g = read_gauges(rundir, model_of(c))
    if g is None:
        return fail("no gauge output")
    a, xw, zb, thr = ck["a"], ck.get("x_wall", 0.0), ck.get("z_bed", 0.0), ck["h_thr"]
    gg = 9.81
    data = sorted(read_refdata(ck["data"]), key=lambda r: r[1])  # (T, Z) sorted by Z
    Zd, Td = [r[1] for r in data], [r[0] for r in data]
    txt = ["| Z = x/a | T sim | T exp | rel. error |", "|---|---|---|---|"]
    errs = []
    for i, x in enumerate(g[0]):
        Z = (x - xw) / a
        if Z < Zd[0] or Z > Zd[-1]:
            continue
        ta = None
        for t, y in zip(g[2], g[3][i]):
            if y - zb > thr and y < 1e10:
                ta = t
                break
        Te = interp(Zd, Td, Z)
        if ta is None:
            txt.append("| %.2f | not reached | %.3f | |" % (Z, Te))
            errs.append(float("inf"))
            continue
        Ts = ta * math.sqrt(2 * gg / a)
        errs.append(abs(Ts - Te) / Te)
        txt.append("| %.2f | %.3f | %.3f | %+.3f |" % (Z, Ts, Te, (Ts - Te) / Te))
    if not errs:
        return fail("no gauge inside the measured range")
    e, em = max(errs), sum(errs) / len(errs)
    ok = e <= ck["tol_max"] and em <= ck["tol_mean"]
    return result(ok, e, "max %g / mean %g" % (ck["tol_max"], ck["tol_mean"]),
                  "\n".join(txt + ["", "mean rel. error %.3f" % em]), error_mean=em)


def check_vortex_shedding(c, rundir):
    """Laminar flow past a circular cylinder: Strouhal number from the lift-force zero crossings and
    mean drag coefficient C_D = F_x / (0.5 rho U^2 D W) in the window, against the reference values
    (Williamson's St(Re) fit; C_D from the compilation of Qu et al. 2013)."""
    ck = c["check"]
    f = read_force(rundir, model_of(c), ck.get("force", 1))
    if f is None or not f["t"]:
        return fail("no force output (P 81)")
    t0, t1 = ck["t_start"], ck["t_end"]
    tt = [t for t in f["t"] if t0 <= t <= t1]
    lift = window(f["t"], f[ck.get("lift", "Fz")], t0, t1)
    drag = window(f["t"], f["Fx"], t0, t1)
    if len(tt) < 10:
        return fail("force output does not cover the window")
    zc = zero_up_crossings(tt, lift)
    U, D, W, rho = ck["U"], ck["D"], ck["W"], ck.get("rho", 1000.0)
    q = 0.5 * rho * U * U * D * W
    Cd = sum(drag) / len(drag) / q
    Cl = (max(lift) - min(lift)) / 2 / q
    if len(zc) < 3 or Cl < ck.get("min_Cl", 0.05):
        r = fail("no periodic shedding in the window (%d lift crossings, C_L amplitude %.2g, C_D %.3f)"
                 % (len(zc), Cl, Cd))
        r.update(Cd=Cd, Cl=Cl)
        return r
    T = (zc[-1] - zc[0]) / (len(zc) - 1)
    St = D / (U * T)
    Re = U * D / ck["nu"]
    A, B, C = -3.3265, 0.1816, 1.6e-4   # Williamson fit (Beaudan & Moin 1994, Eq. 1), 50 < Re < 180
    St_ref = ck.get("St_ref", A / Re + B + C * Re)
    eSt = (St - St_ref) / St_ref
    eCd = (Cd - ck["Cd_ref"]) / ck["Cd_ref"]
    ok = abs(eSt) <= ck["tol_St"] and abs(eCd) <= ck["tol_Cd"]
    txt = ("Re = %.1f, %d shedding periods in %.1f..%.1f s\n\n"
           "| quantity | REEF3D | reference | rel. error | tol |\n|---|---|---|---|---|\n"
           "| St | %.4f | %.4f | %+.3f | ±%g |\n| C_D (mean) | %.3f | %.3f | %+.3f | ±%g |\n"
           "| C_L (amplitude) | %.3f | (rms 0.225-0.235) | | |") % (
        Re, len(zc) - 1, t0, t1, St, St_ref, eSt, ck["tol_St"], Cd, ck["Cd_ref"], eCd, ck["tol_Cd"], Cl)
    return result(ok, max(abs(eSt), abs(eCd)), "St ±%g / Cd ±%g" % (ck["tol_St"], ck["tol_Cd"]), txt,
                  St=St, Cd=Cd, Cl=Cl)


def hulme_prediction(table, mass_ratio, R, g=9.81):
    """Linear heave of a floating hemisphere-shaped body (half-submerged sphere): natural frequency
    from omega^2 (m + a33(omega)) = rho g pi R^2 with Hulme's (1982) added mass, damping ratio from
    b33. A = a33/(rho V), B = b33/(rho V omega), V = 2/3 pi R^3, m = mass_ratio * rho V.
    -> (period, damping ratio zeta, kR)"""
    kR, A, Bd = [r[0] for r in table], [r[1] for r in table], [r[2] for r in table]
    k = 1.0
    for _ in range(100):
        a = interp(kR, A, k)
        k_new = 1.5 / (mass_ratio + a)       # omega^2 R/g = (pi R^2 R) / ((m/rho + a V)) / R ... = 1.5/(mr + A)
        if abs(k_new - k) < 1e-10:
            break
        k = k_new
    a, b = interp(kR, A, k), interp(kR, Bd, k)
    omega = math.sqrt(k * g / R)
    zeta = b / (2.0 * (mass_ratio + a))
    return 2 * math.pi / omega * 1.0 / math.sqrt(1 - zeta * zeta), zeta, k


def check_heave_decay(c, rundir):
    """Free heave decay of a floating sphere of half the water density (CFD 6DOF). References:
    - equilibrium: Archimedes, the centre settles at the still water level (z_eq);
    - damped period and damping ratio from linear potential theory with Hulme's (1982)
      hemisphere added mass and radiation damping (refdata/misc/hulme_hemisphere_heave.txt).
    The period is taken from the first n_periods periods of z(t) - z_eq (successive maxima),
    the damping ratio from the logarithmic decrement of the first peaks."""
    ck = c["check"]
    pos = read_6dof_position(rundir, model_of(c))
    if pos is None or len(pos["t"]) < 10:
        return fail("no 6DOF position output")
    zeq = ck["z_eq"]
    ts, zs = pos["t"], [z - zeq for z in pos["z"]]
    ex = [e for e in extrema(ts, zs) if abs(e[1]) > ck.get("min_amp", 0.0)]
    peaks = [e for e in ex if e[2] == "max"]
    troughs = [e for e in ex if e[2] == "min"]
    allp = sorted(peaks + troughs)
    if len(allp) < 4:
        return fail("less than two oscillations in z(t)")
    n = min(ck.get("n_periods", 2), (len(allp) - 1) // 2)
    T = (allp[2 * n][0] - allp[0][0]) / n
    # log decrement from successive extrema of |z|: one half period apart
    dec = [math.log(abs(allp[i][1]) / abs(allp[i + 1][1])) for i in range(min(len(allp) - 1, 2 * n))
           if allp[i + 1][1] != 0]
    delta_half = sum(dec) / len(dec)
    zeta = (2 * delta_half) / math.sqrt(4 * math.pi ** 2 + (2 * delta_half) ** 2)
    table = read_refdata(ck["hulme"])
    T_ref, zeta_ref, kR = hulme_prediction(table, ck["mass_ratio"], ck["R"])
    tail = [z for t, z in zip(ts, pos["z"]) if t >= ts[-1] - ck.get("eq_window", 1.0)]
    z_end = sum(tail) / len(tail)
    eT = (T - T_ref) / T_ref
    ez = (zeta - zeta_ref) / zeta_ref
    eeq = abs(z_end - zeq) / ck["R"]
    ok = abs(eT) <= ck["tol_T"] and abs(ez) <= ck["tol_zeta"] and eeq <= ck["tol_eq"]
    txt = ("Linear theory (Hulme 1982): kR = %.3f\n\n"
           "| quantity | REEF3D | reference | error | tol |\n|---|---|---|---|---|\n"
           "| damped period T [s] | %.4f | %.4f | %+.3f | ±%g |\n"
           "| damping ratio ζ | %.4f | %.4f | %+.3f | ±%g |\n"
           "| equilibrium z [m] (end of run) | %.4f | %.4f | %.3f R | %g R |\n\n"
           "first extrema of z - z_eq: %s") % (
        kR, T, T_ref, eT, ck["tol_T"], zeta, zeta_ref, ez, ck["tol_zeta"], z_end, zeq, eeq, ck["tol_eq"],
        ", ".join("%.3f s: %+.4f" % (e[0], e[1]) for e in allp[:6]))
    return result(ok, abs(eT), "T ±%g / ζ ±%g / eq %g R" % (ck["tol_T"], ck["tol_zeta"], ck["tol_eq"]), txt,
                  T=T, zeta=zeta)


def check_runup_law(c, rundir):
    """Maximum run-up of a non-breaking solitary wave on a plane beach against Synolakis' (1987)
    run-up law R/d = 2.831 sqrt(cot beta) (H/d)^(5/4). R from the NHFLOW run-up gauge (P 134,
    elevation of the wet-dry front) or from SFLOW wave gauges on the beach (highest wet gauge)."""
    ck = c["check"]
    d, Hd, cot = ck["d"], ck["H_over_d"], ck["cot_beta"]
    R_ref = 2.831 * math.sqrt(cot) * Hd ** 1.25 * d
    model = model_of(c)
    R = None
    fn = os.path.join(rundir, "REEF3D_NHFLOW_RUNUP", "REEF3D-NHFLOW-runup-max-x.dat")
    if model == "NHFLOW" and os.path.exists(fn):
        rows = [r for r in numeric_rows(fn) if len(r) == 3]
        if rows:
            R = rows[-1][2] - ck["swl"]
    else:
        fn = os.path.join(rundir, "REEF3D_%s_RUNUP" % model, "REEF3D-%s-runup-x.dat" % model)
        rows = [r for r in numeric_rows(fn) if len(r) == 3] if os.path.exists(fn) else []
        if rows:
            R = max(r[2] for r in rows) - ck["swl"]
    if R is None:
        return fail("no run-up output")
    e = (R - R_ref) / R_ref
    ok = abs(e) <= ck["tol"]
    return result(ok, abs(e), ck["tol"], "H/d = %g, cot β = %g: R/d REEF3D %.4f, run-up law %.4f, rel. error %+.3f" % (
        Hd, cot, R / d, R_ref / d, e), R=R)


def check_profiles(c, rundir):
    """Free-surface profiles at given non-dimensional times (Synolakis 1987, H/d = 0.28) from the
    WSFLINE output (P 52 with P 55 interval). x and eta are scaled by d, t by sqrt(d/g). The time
    origin of the measured profiles is fitted on the first profile (shift t0 of the simulation).
    Error per profile: rms(eta_sim - eta_exp)/H at the measured points."""
    ck = c["check"]
    lines = read_wsflines(rundir, model_of(c))
    if not lines:
        return fail("no WSFLINE output (P 52)")
    d, H, g = ck["d"], ck["H_over_d"] * ck["d"], 9.81
    tsc = math.sqrt(d / g)
    swl = ck["swl"]
    x0 = ck.get("x_shore", 0.0)  # x of the still-water shoreline in the simulation
    profs = [(tstar, read_refdata(fn)) for tstar, fn in ck["profiles"]]

    def prof_err(line, rows):
        _, xs, es = line
        dd = []
        for xd, ed in rows:
            v = interp(xs, es, x0 + xd * d)
            if v == v:
                dd.append((v - swl) / d - ed)
        return rms(dd) * d / H if len(dd) > 0.7 * len(rows) else float("inf")

    # fit the time origin on the first profile
    t1, r1 = profs[0]
    best = min(lines, key=lambda L: prof_err(L, r1))
    t0 = best[0] - t1 * tsc
    ok, worst = True, 0.0
    txt = ["time origin fitted on t* = %g: t0 = %.3f s" % (t1, t0), "",
           "| t* | sim time [s] | rel. rms error (/H) | |", "|---|---|---|---|"]
    for ts, rows in profs:
        tsim = t0 + ts * tsc
        L = min(lines, key=lambda L: abs(L[0] - tsim))
        e = prof_err(L, rows) if abs(L[0] - tsim) < 0.05 * tsc else float("inf")
        good = e <= ck["tol"]
        ok &= good
        worst = max(worst, e)
        txt.append("| %g | %.3f (line at %.3f) | %.3f | %s |" % (ts, tsim, L[0], e, "" if good else "**FAIL**"))
    return result(ok, worst, ck["tol"], "\n".join(txt), t0=t0)


def check_section_heights(c, rundir):
    """Wave height along measurement sections (Berkhoff et al. 1982 elliptic shoal): wave gauges
    (P 51) along each section, H = mean zero-crossing height in the window, normalised by the
    incident height H0. Error: rms over all points of (H/H0)_sim - (H/H0)_exp, and per section."""
    ck = c["check"]
    g = read_gauges(rundir, model_of(c))
    if g is None:
        return fail("no gauge output")
    t0, t1, H0 = ck["t_start"], ck["t_end"], ck["H0"]
    allerr = []
    txt = ["| section | points | rms error of H/H0 | max error | |", "|---|---|---|---|---|"]
    ok = True
    for sec in ck["sections"]:
        rows = read_refdata(sec["data"])
        errs = []
        for s, hd in rows:
            x = sec["x"] if "x" in sec else sec["sign"] * s + sec.get("offset", 0.0)
            y = sec["y"] if "y" in sec else sec["sign"] * s + sec.get("offset", 0.0)
            i = gauge_index(g, {"x": x + ck.get("dx", 0.0), "y": y + ck.get("dy", 0.0)})
            h = wave_height(g[2], g[3][i], t0, t1) / H0
            errs.append(h - hd)
        e = rms(errs)
        good = e <= ck["tol"]
        ok &= good
        allerr += errs
        txt.append("| %s | %d | %.3f | %.3f | %s |" % (sec["name"], len(errs), e, max(abs(x) for x in errs),
                                                    "" if good else "**FAIL**"))
    return result(ok, rms(allerr), ck["tol"], "\n".join(txt))


def check_breaking_point(c, rundir):
    """Breaking point of regular waves on a slope (Ting & Kirby): dense wave gauges along the
    slope, H(x) = mean zero-crossing height in the window; the breaking point is the location of
    the maximum wave height. Error: |x_b,sim - x_b,exp| (m)."""
    ck = c["check"]
    g = read_gauges(rundir, model_of(c))
    if g is None:
        return fail("no gauge output")
    t0, t1 = ck["t_start"], ck["t_end"]
    Hs = [(x, wave_height(g[2], g[3][i], t0, t1)) for i, x in enumerate(g[0])]
    Hs.sort()
    xb, Hb = max(Hs, key=lambda r: r[1] if r[1] == r[1] else -1)
    xb_ref = ck["x_b"] + ck.get("x_offset", 0.0)
    e = xb - xb_ref
    ok = abs(e) <= ck["tol"]
    txt = ("breaking point (max H): x = %.3f m (H = %.4f m); Ting & Kirby: x = %.3f m -> %+.3f m (tol ±%g m)\n\n"
           "H(x): %s") % (xb, Hb, xb_ref, e, ck["tol"], ", ".join("%.2f:%.3f" % r for r in Hs))
    return result(ok, abs(e), ck["tol"], txt, x_b=xb, H_b=Hb)


def dambreak_exact(x, t, x0, hl, hr, g=9.81):
    """1D dam break on a flat frictionless bed: Ritter (hr = 0) and Stoker (hr > 0), as in
    SWASHES (Delestre et al. 2013). Returns the water depth."""
    if t <= 0:
        return hl if x <= x0 else hr
    cl = math.sqrt(g * hl)
    if hr <= 0:
        xa, xb = x0 - t * cl, x0 + 2 * t * cl
        if x <= xa:
            return hl
        if x <= xb:
            return 4.0 / (9 * g) * (cl - (x - x0) / (2 * t)) ** 2
        return 0.0
    # Stoker: cm from  -8 g hr cm^2 (sqrt(g hl) - cm)^2 + (cm^2 - g hr)^2 (cm^2 + g hr) = 0
    f = lambda cm: -8 * g * hr * cm * cm * (cl - cm) ** 2 + (cm * cm - g * hr) ** 2 * (cm * cm + g * hr)
    lo, hi = math.sqrt(g * hr), cl
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(lo) * f(mid) <= 0:
            hi = mid
        else:
            lo = mid
    cm = 0.5 * (lo + hi)
    xa = x0 - t * cl
    xb = x0 + t * (2 * cl - 3 * cm)
    xc = x0 + t * (2 * cm * cm * (cl - cm)) / (cm * cm - g * hr)
    if x <= xa:
        return hl
    if x <= xb:
        return 4.0 / (9 * g) * (cl - (x - x0) / (2 * t)) ** 2
    if x <= xc:
        return cm * cm / g
    return hr


def check_swe_dambreak(c, rundir):
    """1D dam break (SFLOW) against the exact shallow-water solution (Ritter dry bed, Stoker wet
    bed): depth time series at the wave gauges. Error: sum|h - h_exact| / sum h_exact over all
    gauges and output times in the window (relative L1)."""
    ck = c["check"]
    g = read_gauges(rundir, model_of(c))
    if g is None:
        return fail("no gauge output")
    zb = ck.get("z_bed", 0.0)
    t0, t1 = ck["t_start"], ck["t_end"]
    num = den = 0.0
    txt = ["| x | rel. L1 error |", "|---|---|"]
    for i, x in enumerate(g[0]):
        n_ = d_ = 0.0
        for t, y in zip(g[2], g[3][i]):
            if t0 <= t <= t1:
                h = max(y - zb, 0.0) if y > -1e10 else 0.0
                he = dambreak_exact(x, t, ck["x0"], ck["hl"], ck["hr"])
                n_ += abs(h - he)
                d_ += he
        num += n_
        den += d_
        txt.append("| %.2f | %.4f |" % (x, n_ / d_ if d_ else float("nan")))
    e = num / den
    return result(e <= ck["tol"], e, ck["tol"], "\n".join(txt))


def check_solitary(c, rundir):
    """Solitary wave on constant depth: celerity from the crest arrival times at two gauges against
    c = sqrt(g (d + H)) (Boussinesq/KdV, e.g. Synolakis 1987), and crest height retention H_2/H."""
    ck = c["check"]
    g = read_gauges(rundir, model_of(c))
    if g is None:
        return fail("no gauge output")
    i1, i2 = gauge_index(g, {"x": ck["x1"]}), gauge_index(g, {"x": ck["x2"]})
    e1, e2 = eta_series(c, g, i1), eta_series(c, g, i2)
    k1 = max(range(len(e1)), key=lambda k: e1[k])
    k2 = max(range(len(e2)), key=lambda k: e2[k])
    tc1 = [x for x in extrema(g[2], e1) if x[2] == "max" and abs(x[0] - g[2][k1]) < 0.5]
    tc2 = [x for x in extrema(g[2], e2) if x[2] == "max" and abs(x[0] - g[2][k2]) < 0.5]
    t1_, a1 = (tc1[0][0], tc1[0][1]) if tc1 else (g[2][k1], e1[k1])
    t2_, a2 = (tc2[0][0], tc2[0][1]) if tc2 else (g[2][k2], e2[k2])
    cs = (g[0][i2] - g[0][i1]) / (t2_ - t1_)
    d, H = ck["d"], ck["H"]
    cr = math.sqrt(9.81 * (d + H))
    ec = (cs - cr) / cr
    ea = a2 / H - 1
    ok = abs(ec) <= ck["tol_c"] and abs(ea) <= ck["tol_H"]
    txt = ("| quantity | REEF3D | reference | rel. error | tol |\n|---|---|---|---|---|\n"
           "| celerity [m/s] | %.4f | %.4f | %+.4f | ±%g |\n| crest height at x = %.1f m [m] | %.4f | %.4f | %+.4f | ±%g |\n\n"
           "crest at x = %.1f m: %.4f m") % (cs, cr, ec, ck["tol_c"], g[0][i2], a2, H, ea, ck["tol_H"], g[0][i1], a1)
    return result(ok, max(abs(ec), abs(ea)), "c ±%g / H ±%g" % (ck["tol_c"], ck["tol_H"]), txt)


def check_scour_depth(c, rundir):
    """Equilibrium local scour depth at a vertical circular pier, S/D against the band spanned by
    established design formulas (HEC-18/CSU; Sumer et al. 1992 S/D = 1.3, sigma 0.7). Uses the
    maximum bed change gauge (P 122)."""
    ck = c["check"]
    fn = os.path.join(rundir, "REEF3D_CFD_Sediment", "REEF3D-CFD-Sediment-Max.dat")
    if not os.path.exists(fn):
        return fail("no sediment max output (P 122)")
    rows = [r for r in numeric_rows(fn) if len(r) == 2]   # sediment time, lowest bed level (bedprobe_max)
    if not rows:
        return fail("empty sediment max output")
    S = ck.get("bed0", 0.0) - rows[-1][1]
    SD = S / ck["D"]
    lo, hi = ck["band"]
    ok = lo <= SD <= hi
    return result(ok, SD, "%g <= S/D <= %g" % (lo, hi),
                  "S/D = %.3f at sediment time %.0f s (band %g..%g; HEC-18 CSU estimate %.2f, Sumer et al. 1.3)" % (
                      SD, rows[-1][0], lo, hi, ck.get("hec18", float("nan"))), SD=SD)


CHECKS = {
    "timeseries": check_timeseries,
    "theory_gauges": check_theory_gauges,
    "front_arrival": check_front_arrival,
    "vortex_shedding": check_vortex_shedding,
    "heave_decay": check_heave_decay,
    "runup_law": check_runup_law,
    "profiles": check_profiles,
    "section_heights": check_section_heights,
    "breaking_point": check_breaking_point,
    "swe_dambreak": check_swe_dambreak,
    "solitary": check_solitary,
    "scour_depth": check_scour_depth,
}


# ==============================================================================================
# report and plots
# ==============================================================================================

def evaluate(names, outroot, level):
    results = []
    for n in names:
        c = load_case(n, level)
        rd = os.path.join(outroot, n)
        st = run_status(rd)
        if st.get("status") != "ok":
            results.append((c, fail(st.get("status"))))
            continue
        try:
            r = CHECKS[c["check"]["type"]](c, rd)
        except Exception as ex:  # a broken output file must not stop the report
            r = fail("checker error: %r" % ex)
        r["time_reef3d"] = st.get("time_reef3d")
        if c.get("xfail"):
            r["xfail"] = c["xfail"]
        results.append((c, r))
    lines = ["# REEF3D literature benchmarks (%s level)" % level, "", "run: `%s`" % outroot, "",
             "| case | model | reference | result | error | tolerance | REEF3D time [s] |",
             "|---|---|---|---|---|---|---|"]
    for c, r in results:
        lines.append("| %s | %s | %s | %s | %s | %s | %s |" % (
            c["name"], model_of(c), c.get("reference", {}).get("short", ""),
            verdict(r), "%.3g" % r["error"] if r["error"] == r["error"] else "–",
            r.get("tol", "") if r.get("note") is None else r["note"], r.get("time_reef3d", "")))
    lines += ["", "XFAIL: known deviation from the reference (reason in the details), does not fail the suite; "
              "XPASS: a known deviation passes now - check and remove the xfail entry.", "", "## Details", ""]
    for c, r in results:
        lines += ["### %s" % c["name"], "", c.get("description", ""), ""]
        if r.get("xfail"):
            lines += ["**Known deviation (xfail):** %s" % r["xfail"], ""]
        if c.get("check", {}).get("note"):
            lines += ["Note: %s" % c["check"]["note"], ""]
        lines += [
                  "Reference: %s" % c.get("reference", {}).get("citation", ""), "", r.get("rows_text", ""), ""]
    txt = "\n".join(lines)
    with open(os.path.join(outroot, "benchmark.md"), "w") as f:
        f.write(txt + "\n")
    with open(os.path.join(outroot, "benchmark.json"), "w") as f:
        json.dump([{"case": c["name"], **{k: v for k, v in r.items() if k != "rows_text"}} for c, r in results],
                  f, indent=1, default=str)
    print(txt.split("\n## Details")[0])
    print("\nreport: %s" % os.path.join(outroot, "benchmark.md"))
    return 0 if all(r["ok"] or (r.get("xfail") and r["error"] == r["error"]) for _, r in results) else 1


def verdict(r):
    if r.get("xfail"):
        return "XPASS" if r["ok"] else "XFAIL"
    return "PASS" if r["ok"] else "**FAIL**"


def plot_case(c, rundir, out_png):
    """overview plot of simulated vs reference data (the checker's inputs)"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ck = c["check"]
    typ = ck["type"]
    model = model_of(c)
    if typ == "timeseries":
        res = check_timeseries(c, rundir)
        lag = res.get("lag", 0.0)
        n = len(ck["signals"])
        fig, axs = plt.subplots(n, 1, figsize=(8, 1.8 * n + 0.6), sharex=False, squeeze=False)
        for ax, s in zip(axs[:, 0], ck["signals"]):
            sim = sim_signal(c, rundir, s)
            ref = ref_signal(s, ck.get("window"))
            ax.plot(ref[0], ref[1], "o", ms=3, color="k", label="measured")
            ts = [t - lag for t in sim[0]]
            t0, t1 = ref[0][0], ref[0][-1]
            ax.plot([t for t in ts if t0 - 1 <= t <= t1 + 1], [y for t, y in zip(ts, sim[1]) if t0 - 1 <= t <= t1 + 1],
                    "-", color="C3", lw=1.2, label="REEF3D")
            ax.set_ylabel(s.get("name"))
        axs[0, 0].legend(loc="upper right", fontsize=8)
        axs[-1, 0].set_xlabel("t (data time base) [s]")
    elif typ == "theory_gauges":
        g, th = read_gauges(rundir, model), read_gauges(rundir, model, theory=True)
        n = len(g[0])
        fig, axs = plt.subplots(n, 1, figsize=(8, 1.6 * n + 0.6), squeeze=False)
        for i, ax in enumerate(axs[:, 0]):
            ax.plot(th[2], th[3][i], "-", color="k", lw=1, label="theory")
            ax.plot(g[2], eta_series(c, g, i), "-", color="C3", lw=1, label="REEF3D")
            ax.set_xlim(ck["t_start"], ck["t_end"])
            ax.set_ylabel("x=%.1f" % g[0][i])
        axs[0, 0].legend(fontsize=8)
    elif typ in ("front_arrival",):
        g = read_gauges(rundir, model)
        a = ck["a"]
        data = read_refdata(ck["data"])
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot([r[0] for r in data], [r[1] for r in data], "o", color="k", label="Martin & Moyce")
        T, Z = [], []
        for i, x in enumerate(g[0]):
            for t, y in zip(g[2], g[3][i]):
                if y - ck.get("z_bed", 0.0) > ck["h_thr"] and y < 1e10:
                    T.append(t * math.sqrt(2 * 9.81 / a))
                    Z.append((x - ck.get("x_wall", 0.0)) / a)
                    break
        ax.plot(T, Z, "s-", color="C3", label="REEF3D")
        ax.set_xlabel("T = t sqrt(2g/a)")
        ax.set_ylabel("Z = x/a")
        ax.legend()
    elif typ == "heave_decay":
        pos = read_6dof_position(rundir, model)
        fig, ax = plt.subplots(figsize=(8, 3.5))
        ax.plot(pos["t"], pos["z"], color="C3", label="REEF3D")
        ax.axhline(ck["z_eq"], color="k", lw=0.8, ls="--", label="Archimedes equilibrium")
        ax.set_xlabel("t [s]")
        ax.set_ylabel("z_G [m]")
        ax.legend()
    elif typ in ("section_heights", "breaking_point", "swe_dambreak", "solitary", "vortex_shedding", "runup_law",
                 "profiles", "scour_depth"):
        fig, ax = plt.subplots(figsize=(8, 4))
        if typ == "section_heights":
            g = read_gauges(rundir, model)
            for k, sec in enumerate(ck["sections"]):
                rows = read_refdata(sec["data"])
                ax.plot([r[0] + 10 * k for r in rows], [r[1] for r in rows], "o", ms=3, color="k")
                hs = []
                for s, _ in rows:
                    x = sec["x"] if "x" in sec else sec["sign"] * s + sec.get("offset", 0.0)
                    y = sec["y"] if "y" in sec else sec["sign"] * s + sec.get("offset", 0.0)
                    i = gauge_index(g, {"x": x + ck.get("dx", 0.0), "y": y + ck.get("dy", 0.0)})
                    hs.append(wave_height(g[2], g[3][i], ck["t_start"], ck["t_end"]) / ck["H0"])
                ax.plot([r[0] + 10 * k for r in rows], hs, "-", color="C3")
                ax.text(10 * k - 2, 2.2, sec["name"])
            ax.set_ylabel("H/H0 (sections offset by 10 m)")
        elif typ == "breaking_point":
            g = read_gauges(rundir, model)
            Hs = sorted((x, wave_height(g[2], g[3][i], ck["t_start"], ck["t_end"])) for i, x in enumerate(g[0]))
            ax.plot([r[0] for r in Hs], [r[1] for r in Hs], "-", color="C3", label="REEF3D H(x)")
            ax.axvline(ck["x_b"] + ck.get("x_offset", 0.0), color="k", ls="--", label="Ting & Kirby x_b")
            ax.set_xlabel("x [m]")
            ax.legend()
        elif typ == "swe_dambreak":
            g = read_gauges(rundir, model)
            for i, x in enumerate(g[0]):
                ax.plot(g[2], [max(y - ck.get("z_bed", 0.0), 0) for y in g[3][i]], color="C%d" % i, lw=1)
                ax.plot(g[2], [dambreak_exact(x, t, ck["x0"], ck["hl"], ck["hr"]) for t in g[2]], "--", color="C%d" % i,
                        lw=1)
            ax.set_xlabel("t [s]")
            ax.set_ylabel("h [m] (solid REEF3D, dashed exact)")
        elif typ == "solitary":
            g = read_gauges(rundir, model)
            for i, x in enumerate(g[0]):
                ax.plot(g[2], eta_series(c, g, i), label="x=%.1f" % x)
            ax.legend()
        elif typ == "vortex_shedding":
            f = read_force(rundir, model, ck.get("force", 1))
            ax.plot(f["t"], f["Fx"], label="Fx")
            ax.plot(f["t"], f[ck.get("lift", "Fz")], label="lift")
            ax.legend()
        elif typ == "profiles":
            lines = read_wsflines(rundir, model)
            res = check_profiles(c, rundir)
            d = ck["d"]
            tsc = math.sqrt(d / 9.81)
            for k, (ts, fn) in enumerate(ck["profiles"]):
                rows = read_refdata(fn)
                L = min(lines, key=lambda L: abs(L[0] - (res["t0"] + ts * tsc)))
                ax.plot([r[0] for r in rows], [r[1] + 0.3 * k for r in rows], "o", ms=2, color="k")
                ax.plot([(x - ck.get("x_shore", 0.0)) / d for x in L[1]], [(e - ck["swl"]) / d + 0.3 * k for e in L[2]],
                        color="C3", lw=1)
            ax.set_xlim(-20, 5)
            ax.set_xlabel("x/d")
            ax.set_ylabel("eta/d (offset 0.3 per profile)")
        else:
            ax.text(0.1, 0.5, "see benchmark.md", transform=ax.transAxes)
    else:
        return False
    fig.suptitle(c["name"], fontsize=10)
    fig.tight_layout()
    fig.savefig(out_png, dpi=110)
    plt.close(fig)
    return True


# ==============================================================================================

def suite_level(outroot, default):
    try:
        with open(os.path.join(outroot, "suite.json")) as f:
            return json.load(f).get("level", default)
    except OSError:
        return default


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def sel(p, level=True):
        p.add_argument("--cases", nargs="*")
        p.add_argument("--tags", nargs="*")
        if level:
            p.add_argument("--level", choices=LEVELS, default="nightly")

    p = sub.add_parser("list"); sel(p)
    p = sub.add_parser("write"); sel(p); p.add_argument("--out", required=True)
    p.add_argument("--keep", action="store_true", help="keep the VTU/VTP print keys")
    p = sub.add_parser("run"); sel(p)
    p.add_argument("--reef3d", required=True); p.add_argument("--divemesh", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--mpirun", default=os.environ.get("REEF3D_MPIRUN", "mpirun"))
    p.add_argument("--timeout", type=int, default=6 * 3600)
    p.add_argument("--keep", action="store_true", help="keep grids and VTU/VTP output")
    p.add_argument("--max-np", type=int, help="cap the number of MPI ranks (results change slightly with the "
                                             "decomposition)")
    p.add_argument("--steps", type=int, help="smoke test: stop after this many time steps (N 45)")
    p.add_argument("--check", action="store_true")
    p = sub.add_parser("check"); sel(p, False); p.add_argument("run")
    p = sub.add_parser("plot"); sel(p, False); p.add_argument("run")
    a = ap.parse_args()

    if a.cmd == "list":
        for n in select_cases(a.cases, a.tags, a.level):
            c = load_case(n, a.level)
            print("%-38s %-6s np=%-3d t=%-7s %s\n    %s" % (n, model_of(c), c["np"], c.get("time"),
                                                         ",".join(c["tags"]), c.get("reference", {}).get("short", "")))
        return 0
    if a.cmd == "write":
        for n in select_cases(a.cases, a.tags, a.level):
            write_inputs(load_case(n, a.level), os.path.join(a.out, n), a.keep)
            print(os.path.join(a.out, n))
        return 0
    if a.cmd == "run":
        names = select_cases(a.cases, a.tags, a.level)
        if not names:
            sys.exit("no cases selected")
        os.makedirs(a.out, exist_ok=True)
        with open(os.path.join(a.out, "suite.json"), "w") as f:
            json.dump({"level": a.level, "cases": names, "reef3d": os.path.abspath(a.reef3d),
                       "date": time.strftime("%Y-%m-%d %H:%M:%S"), "host": os.uname().nodename}, f, indent=1)
        failed = 0
        for n in names:
            c = load_case(n, a.level)
            if a.max_np:
                c["np"] = min(c["np"], a.max_np)
            if a.steps:  # smoke test: a few time steps only (the checks will not be meaningful)
                c["ctrl_lines"] = apply_overrides(c["ctrl_lines"], {"N 45": str(a.steps)}, [])
            print("%-38s np=%-3d ... " % (n, c["np"]), end="", flush=True)
            info = run_case(c, a, a.out)
            print("%s (%.0f s)" % (info["status"], info["time_divemesh"] + info["time_reef3d"]), flush=True)
            failed += info["status"] != "ok"
        rc = 1 if failed else 0
        if a.check:
            rc |= evaluate(names, a.out, a.level)
        return rc
    if a.cmd in ("check", "plot"):
        level = suite_level(a.run, "nightly")
        names = [n for n in select_cases(a.cases, a.tags, level) if os.path.isdir(os.path.join(a.run, n))]
        if not names:
            sys.exit("no case folders of the selected cases in %s" % a.run)
        if a.cmd == "check":
            return evaluate(names, a.run, level)
        for n in names:
            c = load_case(n, level)
            rd = os.path.join(a.run, n)
            if run_status(rd).get("status") != "ok":
                continue
            png = os.path.join(rd, n + ".png")
            try:
                if plot_case(c, rd, png):
                    print(png)
            except Exception as ex:
                print("%s: plot failed (%r)" % (n, ex))
        return 0


if __name__ == "__main__":
    sys.exit(main())
