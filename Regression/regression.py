#!/usr/bin/env python3
"""
REEF3D regression test suite.

Pure Python 3 standard library (no numpy needed).

Typical use
-----------
  # 1. run the suite with a reference binary (e.g. hans_dev before a change)
  ./regression.py run --reef3d /path/REEF3D_ref --divemesh /path/DiveMESH --out runs/ref

  # 2. run the same cases with the modified binary
  ./regression.py run --reef3d /path/REEF3D_new --divemesh /path/DiveMESH --out runs/new

  # 3. compare: bitwise for refactors, tolerance for physics changes
  ./regression.py compare runs/ref runs/new            # requires bitwise identity
  ./regression.py compare runs/ref runs/new --require close --rtol 1e-6

  # or 1-3 in one go
  ./regression.py ab --ref-bin REEF3D_ref --new-bin REEF3D_new --divemesh DiveMESH --out runs/ab

  # long-term: store compact references in cases/<case>/reference.json and check against them
  ./regression.py bless runs/ref
  ./regression.py check runs/new --rtol 1e-6

Case selection: --cases <glob> [<glob> ...] and/or --tags <tag> [...] (e.g. --tags quick).

How it works
------------
Each case lives in cases/<name>/ with a case.json, a DIVEMesh control.txt and a REEF3D ctrl.txt
(or "base": "<other case>" plus key overrides). The runner writes the input files, runs DIVEMesh
and REEF3D with REEF3D_REGRESSION_DIR set, which makes REEF3D (src/regression_dump.cpp) write
  - steps_r<rank>.txt   per-step count, simtime, dt and field norms in hexfloat (exact)
  - state_<count>_r<rank>.bin   full double-precision state (initial and final, optional every n)
The comparison works on these exact dumps and, with a tolerance, on the normal text output
(wave gauges, probes, forces, 6DOF) in the REEF3D_CFD_* folders.
"""

import argparse
import fnmatch
import glob
import json
import math
import os
import re
import shlex
import shutil
import struct
import subprocess
import sys
import time
from array import array

HERE = os.path.dirname(os.path.abspath(__file__))
CASES_DIR = os.path.join(HERE, "cases")
DUMP_DIR = "regression_dump"

# text output folders that are compared numerically (with tolerance); VTU etc. are ignored
TEXT_OUTPUT_GLOBS = [
    "REEF3D_CFD_WSF/*.dat",
    "REEF3D_CFD_ProbePoint/*.dat",
    "REEF3D_CFD_PressureProbe/*.dat",
    "REEF3D_CFD_ProbeLine/*.dat",
    "REEF3D_CFD_Force/*.dat",
    "REEF3D_CFD_6DOF/*.dat",
    "REEF3D_CFD_WSFLINE/*.dat",
    "REEF3D_CFD_CPM_Particle/*.dat",
    "REEF3D_CFD_Acoustics/*.dat",
    "REEF3D_FEM/*.dat",
    "REEF3D_NHFLOW_WSF/*.dat",
    "REEF3D_NHFLOW_ProbePoint/*.dat",
    "REEF3D_NHFLOW_Force/*.dat",
    "REEF3D_NHFLOW_WSFLINE/*.dat",
    "REEF3D_NHFLOW_6DOF/*.dat",
    "REEF3D_NHFLOW_AMR/*.dat",
    "REEF3D_NHFLOW_Particles/*.dat",
    "REEF3D_NHFLOW_Boom/*.dat",
    "REEF3D_DEM/*.dat",
    "REEF3D_SFLOW_AMR/*.dat",
    "REEF3D_FNPF_AMR/*.dat",
    "REEF3D_FNPF_WSF/*.dat",
    "REEF3D_FNPF_ProbePoint/*.dat",
    "REEF3D_FNPF_WSFLINE/*.dat",
    "REEF3D_FNPF_6DOF/*.dat",
    "REEF3D_SFLOW_WSF/*.dat",
    "REEF3D_SFLOW_ProbePoint/*.dat",
    "REEF3D_SFLOW_WSFLINE/*.dat",
    "REEF3D_SFLOW_6DOF/*.dat",
    "REEF3D_SEASTATE_Log/*.dat",
    "REEF3D_SEASTATE_Spectra/*.dat",
]

# output folders deleted after a run unless --keep (large and not compared)
BULK_OUTPUT = ["REEF3D_CFD_VTU", "REEF3D_CFD_6DOF_VTP", "REEF3D_CFD_6DOF_Normals_VTP",
               "REEF3D_CFD_FSF", "DIVEMesh_Paraview",
               "REEF3D_NHFLOW_VTU", "REEF3D_NHFLOW_VTP_FSF", "REEF3D_NHFLOW_VTP_BED",
               "REEF3D_FNPF_VTU", "REEF3D_FNPF_VTP_FSF", "REEF3D_FNPF_VTP_BED",
               "REEF3D_SFLOW_VTP_FSF", "REEF3D_SFLOW_VTP_BED",
               "REEF3D_NHFLOW_6DOF_VTP", "REEF3D_NHFLOW_6DOF_Normals_VTP", "REEF3D_NHFLOW_6DOF_STL",
               "REEF3D_FNPF_6DOF_VTP", "REEF3D_FNPF_6DOF_Normals_VTP", "REEF3D_FNPF_6DOF_STL",
               "REEF3D_CFD_6DOF_STL", "REEF3D_SFLOW_6DOF_VTP", "REEF3D_SEASTATE_VTP",
               "REEF3D_DEM_VTP"]

LEVELS = ["identical", "close", "different", "failed"]


# ----------------------------------------------------------------------------------------------
# case definitions
# ----------------------------------------------------------------------------------------------

# cap on the number of ranks of every case (--max-np): both binaries of an A/B run use the same
# reduced decomposition, so the comparison stays valid on a small machine; the stored references
# (bless/check) are for the case's own np
MAX_NP = 0


def load_case(name, _seen=None):
    """Load a case, resolving 'base' inheritance. Returns dict with resolved control/ctrl lines."""
    _seen = _seen or set()
    if name in _seen:
        raise ValueError("circular base in case %s" % name)
    _seen.add(name)
    cdir = os.path.join(CASES_DIR, name)
    with open(os.path.join(cdir, "case.json")) as f:
        c = json.load(f)
    c["name"] = name
    c["dir"] = cdir
    if c.get("base"):
        b = load_case(c["base"], _seen)
        control = list(b["control_lines"])
        ctrl = list(b["ctrl_lines"])
        files = list(b["files_abs"])
        # inherit everything that is not specific to the base case's own files/overrides
        for k, v in b.items():
            if k not in ("name", "dir", "base", "control_lines", "ctrl_lines", "files_abs", "files",
                         "control_set", "ctrl_set", "control_add", "ctrl_add", "description", "covers"):
                c.setdefault(k, v)
    else:
        control = read_lines(os.path.join(cdir, "control.txt"))
        ctrl = read_lines(os.path.join(cdir, "ctrl.txt"))
        files = []
    control = apply_overrides(control, c.get("control_set", {}), c.get("control_add", []))
    ctrl = apply_overrides(ctrl, c.get("ctrl_set", {}), c.get("ctrl_add", []))
    files += [os.path.join(cdir, x) for x in c.get("files", [])]
    c["control_lines"] = control
    c["ctrl_lines"] = ctrl
    c["files_abs"] = files
    c.setdefault("np", 1)
    if MAX_NP > 0:
        c["np"] = min(c["np"], MAX_NP)
    c.setdefault("steps", 50)
    c.setdefault("dump_every", 0)
    c.setdefault("solver", "cfd")
    c.setdefault("tags", [])
    return c


def read_lines(path):
    with open(path) as f:
        return [l.rstrip("\n") for l in f]


def key_of(line):
    t = line.split()
    if len(t) >= 2 and len(t[0]) == 1 and t[0].isalpha() and t[1].isdigit():
        return t[0] + " " + t[1]
    return None


def apply_overrides(lines, setmap, addlist):
    """setmap: {"N 40": "13"} replaces all lines with that key; value null removes them;
    a list value writes several lines. addlist: lines appended as they are."""
    out = [l for l in lines if key_of(l) not in setmap]
    for k, v in setmap.items():
        if v is None:
            continue
        vals = v if isinstance(v, list) else [v]
        for x in vals:
            out.append("%s %s" % (k, x))
    out += list(addlist)
    return out


def all_case_names():
    return sorted(d for d in os.listdir(CASES_DIR)
                  if os.path.isfile(os.path.join(CASES_DIR, d, "case.json")))


def select_cases(patterns, tags):
    names = all_case_names()
    if patterns:
        names = [n for n in names if any(fnmatch.fnmatch(n, p) for p in patterns)]
    if tags:
        names = [n for n in names if set(tags) & set(load_case(n)["tags"])]
    return names


# ----------------------------------------------------------------------------------------------
# running
# ----------------------------------------------------------------------------------------------

def write_inputs(c, rundir, steps_override=None):
    os.makedirs(rundir, exist_ok=True)
    np_ = c["np"]
    steps = steps_override or c["steps"]
    control = apply_overrides(c["control_lines"], {"M 10": str(np_)}, [])
    # the suite controls length and output (no VTU/state files)
    # length: fixed number of steps, or a simulated time ("time" in case.json, validation cases)
    if c.get("time") and not steps_override:
        nsteps, tmax = "100000000", repr(float(c["time"]))
    else:
        nsteps, tmax = str(steps), "1.0e9"
    suite = {"M 10": str(np_), "N 45": nsteps, "N 41": tmax,
             "P 20": None, "P 30": None, "P 40": None, "P 41": None,
             "P 42": None, "P 12": "10"}
    # "keep_keys": print keys a case needs (e.g. P 40 state files for a chain stage)
    for k in c.get("keep_keys", []):
        suite.pop(k, None)
    ctrl = apply_overrides(c["ctrl_lines"], suite, [])
    with open(os.path.join(rundir, "control.txt"), "w") as f:
        f.write("\n".join(control) + "\n")
    with open(os.path.join(rundir, "ctrl.txt"), "w") as f:
        f.write("\n".join(ctrl) + "\n")
    for src in c["files_abs"]:
        shutil.copy(src, rundir)


def run_case(c, args, outroot):
    rundir = os.path.join(outroot, c["name"])
    if os.path.isdir(rundir):
        shutil.rmtree(rundir)
    write_inputs(c, rundir, args.steps)
    # "chain": earlier runs whose output folders this case reads (hydrodynamic coupling):
    # [{"case": "<case>", "copy": ["REEF3D_FNPF_STATE", ...]}], run in <rundir>/_chain
    for stage in c.get("chain", []):
        sc = load_case(stage["case"])
        sinfo = run_case(sc, args, os.path.join(rundir, "_chain"))
        if sinfo["status"] != "ok":
            info = {"case": c["name"], "np": c["np"], "steps": args.steps or c["steps"],
                    "status": "chain stage %s: %s" % (stage["case"], sinfo["status"])}
            return finish_run(rundir, info, args)
        for d in stage.get("copy", []):
            shutil.copytree(os.path.join(rundir, "_chain", stage["case"], d), os.path.join(rundir, d))
    info = {"case": c["name"], "np": c["np"], "steps": args.steps or c["steps"],
            "reef3d": os.path.abspath(args.reef3d), "divemesh": os.path.abspath(args.divemesh),
            "status": "ok", "time_divemesh": 0.0, "time_reef3d": 0.0}
    env = dict(os.environ)
    env["REEF3D_REGRESSION_DIR"] = DUMP_DIR
    env["REEF3D_REGRESSION_EVERY"] = str(c.get("dump_every", 0))
    env.setdefault("OMP_NUM_THREADS", "1")

    t0 = time.time()
    with open(os.path.join(rundir, "divemesh.log"), "w") as log:
        r = subprocess.run([os.path.abspath(args.divemesh)], cwd=rundir, stdout=log,
                           stderr=subprocess.STDOUT, timeout=args.timeout)
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

    if info["status"] == "ok":
        st = read_steps(rundir)
        if not st or not st.get(0):
            info["status"] = "no regression dump (binary without src/regression_dump?)"
        else:
            info["steps_done"] = len(st[0])
            last = st[0][-1]
            info["final_simtime"] = float.fromhex(last[1])
            # NaN check on the per-step norms
            if any(not math.isfinite(float.fromhex(x)) for x in last[3:]):
                info["status"] = "NaN/Inf in solution"
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


def cmd_run(args):
    names = select_cases(args.cases, args.tags)
    if not names:
        sys.exit("no cases selected")
    os.makedirs(args.out, exist_ok=True)
    with open(os.path.join(args.out, "suite.json"), "w") as f:
        json.dump({"reef3d": os.path.abspath(args.reef3d), "cases": names,
                   "date": time.strftime("%Y-%m-%d %H:%M:%S"),
                   "host": os.uname().nodename}, f, indent=1)
    failed = 0
    t_all = time.time()
    for n in names:
        c = load_case(n)
        print("%-40s np=%d steps=%-4d ... " % (n, c["np"], args.steps or c["steps"]), end="", flush=True)
        info = run_case(c, args, args.out)
        print("%s  (%.1f s)" % (info["status"], info["time_divemesh"] + info["time_reef3d"]))
        if info["status"] != "ok":
            failed += 1
    print("%d cases, %d failed, %.0f s" % (len(names), failed, time.time() - t_all))
    return 1 if failed else 0


# ----------------------------------------------------------------------------------------------
# reading dumps
# ----------------------------------------------------------------------------------------------

def read_steps(rundir):
    """{rank: [[count, simtime, dt, norms...] as hex strings]}"""
    out = {}
    for fn in glob.glob(os.path.join(rundir, DUMP_DIR, "steps_r*.txt")):
        rank = int(re.search(r"steps_r(\d+)\.txt$", fn).group(1))
        rows = []
        with open(fn) as f:
            for l in f:
                if l.startswith("#") or not l.strip():
                    continue
                rows.append(l.split())
        out[rank] = rows
    return out


def state_files(rundir):
    """{count: {rank: path}}"""
    out = {}
    for fn in glob.glob(os.path.join(rundir, DUMP_DIR, "state_*_r*.bin")):
        m = re.search(r"state_(\d+)_r(\d+)\.bin$", fn)
        out.setdefault(int(m.group(1)), {})[int(m.group(2))] = fn
    return out


def read_state(path):
    with open(path, "rb") as f:
        buf = f.read()
    if buf[:8] != b"R3DREG01":
        raise ValueError("bad regression state file %s" % path)
    pos = 8
    rank, size, count = struct.unpack_from("<iii", buf, pos); pos += 12
    simtime, = struct.unpack_from("<d", buf, pos); pos += 8
    nf, = struct.unpack_from("<i", buf, pos); pos += 4
    fields = {}
    order = []
    for _ in range(nf):
        name = buf[pos:pos + 16].split(b"\0")[0].decode(); pos += 16
        n, = struct.unpack_from("<q", buf, pos); pos += 8
        a = array("d")
        a.frombytes(buf[pos:pos + 8 * n]); pos += 8 * n
        if sys.byteorder != "little":
            a.byteswap()
        fields[name] = a
        order.append(name)
    return {"rank": rank, "size": size, "count": count, "simtime": simtime,
            "fields": fields, "order": order}


def merged_final_state(rundir):
    sf = state_files(rundir)
    if not sf:
        return None
    count = max(sf)
    ranks = sorted(sf[count])
    states = [read_state(sf[count][r]) for r in ranks]
    return count, states


# ----------------------------------------------------------------------------------------------
# comparison
# ----------------------------------------------------------------------------------------------

def field_diff(a, b):
    """max abs diff, max |a| and number of differing entries (exact)."""
    maxd = 0.0
    maxa = 0.0
    ndiff = 0
    for x, y in zip(a, b):
        ax = abs(x)
        if ax > maxa:
            maxa = ax
        if x != y:
            if x != x and y != y:  # both NaN
                continue
            ndiff += 1
            d = abs(x - y)
            if not (d <= maxd):  # also catches NaN
                maxd = d if d == d else float("inf")
    return maxd, maxa, ndiff


def compare_states(ref_dir, new_dir, rtol, atol):
    res = {"level": "identical", "fields": {}, "notes": []}
    r = merged_final_state(ref_dir)
    n = merged_final_state(new_dir)
    if r is None or n is None:
        res["level"] = "failed"
        res["notes"].append("missing state dump")
        return res
    if r[0] != n[0]:
        res["level"] = "different"
        res["notes"].append("final step differs: %d vs %d" % (r[0], n[0]))
        return res
    if len(r[1]) != len(n[1]):
        res["level"] = "failed"
        res["notes"].append("different number of ranks")
        return res
    for sr, sn in zip(r[1], n[1]):
        for name in sr["order"]:
            fr = sr["fields"][name]
            fnw = sn["fields"].get(name)
            agg = res["fields"].setdefault(name, {"maxdiff": 0.0, "maxref": 0.0, "ndiff": 0,
                                                  "n": 0, "size_mismatch": False})
            if fnw is None or len(fnw) != len(fr):
                agg["size_mismatch"] = True
                continue
            d, m, nd = field_diff(fr, fnw)
            agg["maxdiff"] = max(agg["maxdiff"], d)
            agg["maxref"] = max(agg["maxref"], m)
            agg["ndiff"] += nd
            agg["n"] += len(fr)
    level = 0
    for name, agg in res["fields"].items():
        if agg["size_mismatch"]:
            lv = 2
            res["notes"].append("%s: size mismatch (flags differ?)" % name)
        elif agg["ndiff"] == 0:
            lv = 0
        elif agg["maxdiff"] <= atol + rtol * agg["maxref"]:
            lv = 1
        else:
            lv = 2
        agg["level"] = LEVELS[lv]
        level = max(level, lv)
    res["level"] = LEVELS[level]
    return res


def compare_steps(ref_dir, new_dir, rtol):
    """first step where the exact per-step record differs (rank 0..n)."""
    r = read_steps(ref_dir)
    n = read_steps(new_dir)
    first = None
    dtmax = 0.0
    for rank in sorted(r):
        rr, nn = r[rank], n.get(rank, [])
        for a, b in zip(rr, nn):
            if a != b:
                if first is None or int(a[0]) < first:
                    first = int(a[0])
                break
        for a, b in zip(rr, nn):
            da, db = float.fromhex(a[2]), float.fromhex(b[2])
            dtmax = max(dtmax, abs(da - db) / max(abs(da), 1e-300))
        if len(rr) != len(nn):
            first = first if first is not None else min(len(rr), len(nn))
    return {"first_diff_step": first, "dt_maxrel": dtmax,
            "steps_ref": len(r.get(0, [])), "steps_new": len(n.get(0, []))}


NUM = re.compile(r"^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$")


def read_numeric_table(path):
    rows = []
    with open(path, errors="replace") as f:
        for l in f:
            t = l.split()
            vals = [float(x) for x in t if NUM.match(x)]
            if vals and len(vals) == len(t):
                rows.append(vals)
    return rows


def compare_text(ref_dir, new_dir, rtol, atol):
    res = {"level": "identical", "files": {}}
    level = 0
    for g in TEXT_OUTPUT_GLOBS:
        for fr in sorted(glob.glob(os.path.join(ref_dir, g))):
            rel = os.path.relpath(fr, ref_dir)
            fn = os.path.join(new_dir, rel)
            if not os.path.exists(fn):
                res["files"][rel] = {"level": "failed", "note": "missing"}
                level = max(level, 3)
                continue
            with open(fr, "rb") as a, open(fn, "rb") as b:
                if a.read() == b.read():
                    continue
            ta, tb = read_numeric_table(fr), read_numeric_table(fn)
            if len(ta) != len(tb):
                res["files"][rel] = {"level": "different", "note": "%d vs %d rows" % (len(ta), len(tb))}
                level = max(level, 2)
                continue
            maxd, maxa = 0.0, 0.0
            for ra, rb in zip(ta, tb):
                for x, y in zip(ra, rb):
                    maxa = max(maxa, abs(x))
                    maxd = max(maxd, abs(x - y))
            lv = 1 if maxd <= atol + rtol * maxa else 2
            res["files"][rel] = {"level": LEVELS[lv], "maxdiff": maxd, "maxref": maxa}
            level = max(level, lv)
    res["level"] = LEVELS[level]
    return res


def run_status(d):
    try:
        with open(os.path.join(d, "run.json")) as f:
            return json.load(f)
    except OSError:
        return {"status": "not run"}


def compare_case(name, ref_root, new_root, rtol, atol):
    rd, nd = os.path.join(ref_root, name), os.path.join(new_root, name)
    out = {"case": name}
    sr, sn = run_status(rd), run_status(nd)
    out["time_ref"] = sr.get("time_reef3d")
    out["time_new"] = sn.get("time_reef3d")
    if sr.get("status") != "ok" or sn.get("status") != "ok":
        out["level"] = "failed"
        out["note"] = "ref: %s / new: %s" % (sr.get("status"), sn.get("status"))
        return out
    out["state"] = compare_states(rd, nd, rtol, atol)
    out["steps"] = compare_steps(rd, nd, rtol)
    out["text"] = compare_text(rd, nd, rtol, atol)
    lv = max(LEVELS.index(out["state"]["level"]), LEVELS.index(out["text"]["level"]))
    if out["steps"]["first_diff_step"] is not None:
        lv = max(lv, 1)
    out["level"] = LEVELS[lv]
    return out


def fmt_e(x):
    return "0" if x == 0 else "%.2e" % x


def report(results, ref_root, new_root, rtol, atol):
    lines = []
    lines.append("# REEF3D regression comparison\n")
    lines.append("ref: `%s`  \nnew: `%s`  \nrtol=%g atol=%g\n" % (ref_root, new_root, rtol, atol))
    lines.append("| case | result | first diff step | worst field (max abs diff / max ref) | text output | time ref/new [s] |")
    lines.append("|---|---|---|---|---|---|")
    for r in results:
        if r["level"] == "failed" and "state" not in r:
            lines.append("| %s | **FAILED** | | %s | | |" % (r["case"], r.get("note", "")))
            continue
        worst = ""
        wd = -1
        for name, agg in r["state"]["fields"].items():
            if agg.get("ndiff", 0) and agg["maxdiff"] / max(agg["maxref"], 1e-300) > wd:
                wd = agg["maxdiff"] / max(agg["maxref"], 1e-300)
                worst = "%s %s / %s (%d cells)" % (name, fmt_e(agg["maxdiff"]), fmt_e(agg["maxref"]), agg["ndiff"])
        worst = worst or "; ".join(r["state"]["notes"]) or "–"
        fd = r["steps"]["first_diff_step"]
        tl = r["text"]["level"]
        lines.append("| %s | %s | %s | %s | %s | %s / %s |" % (
            r["case"], r["level"].upper() if r["level"] != "identical" else "identical",
            "–" if fd is None else fd, worst, tl, r.get("time_ref"), r.get("time_new")))
    # details
    lines.append("\n## Details\n")
    for r in results:
        if r["level"] in ("identical",) or "state" not in r:
            continue
        lines.append("### %s (%s)\n" % (r["case"], r["level"]))
        lines.append("| field | level | max abs diff | max ref | differing entries |")
        lines.append("|---|---|---|---|---|")
        for name, agg in r["state"]["fields"].items():
            lines.append("| %s | %s | %s | %s | %d / %d |" % (name, agg.get("level"), fmt_e(agg["maxdiff"]),
                                                              fmt_e(agg["maxref"]), agg["ndiff"], agg["n"]))
        for f, v in r["text"]["files"].items():
            lines.append("- `%s`: %s %s" % (f, v["level"], v.get("note", "max diff %s (max %s)" % (
                fmt_e(v.get("maxdiff", 0)), fmt_e(v.get("maxref", 0))))))
        lines.append("")
    return "\n".join(lines)


def cmd_compare(args):
    names = args.cases_list or sorted(
        set(os.listdir(args.ref)) & set(os.listdir(args.new)) & set(all_case_names()))
    if args.cases:
        names = [n for n in names if any(fnmatch.fnmatch(n, p) for p in args.cases)]
    results = [compare_case(n, args.ref, args.new, args.rtol, args.atol) for n in names]
    txt = report(results, args.ref, args.new, args.rtol, args.atol)
    out = args.report or os.path.join(args.new, "compare.md")
    with open(out, "w") as f:
        f.write(txt + "\n")
    with open(os.path.splitext(out)[0] + ".json", "w") as f:
        json.dump(results, f, indent=1)
    print(txt.split("\n## Details")[0])
    print("\nreport: %s" % out)
    allowed = LEVELS.index(args.require)
    bad = [r["case"] for r in results if LEVELS.index(r["level"]) > allowed]
    if bad:
        print("NOT %s: %s" % (args.require.upper(), ", ".join(bad)))
        return 1
    print("all %d cases %s" % (len(results), args.require if args.require != "identical" else "bitwise identical"))
    return 0


def cmd_ab(args):
    a = argparse.Namespace(**vars(args))
    rc = 0
    for tag, binary in (("ref", args.ref_bin), ("new", args.new_bin)):
        print("== %s: %s" % (tag, binary))
        a.reef3d = binary
        a.out = os.path.join(args.out, tag)
        rc |= cmd_run(a)
    c = argparse.Namespace(ref=os.path.join(args.out, "ref"), new=os.path.join(args.out, "new"),
                           rtol=args.rtol, atol=args.atol, require=args.require, report=None,
                           cases=args.cases, cases_list=select_cases(args.cases, args.tags))
    return cmd_compare(c) | rc


# ----------------------------------------------------------------------------------------------
# stored references (compact, machine independent within tolerance)
# ----------------------------------------------------------------------------------------------

def summarize(rundir):
    st = merged_final_state(rundir)
    steps = read_steps(rundir)
    summ = {"final_count": st[0], "fields": {}, "text": {}}
    summ["final_simtime"] = float.fromhex(steps[0][-1][1])
    summ["dt"] = [float.fromhex(r[2]) for r in steps[0]]
    for s in st[1]:
        for name in s["order"]:
            a = s["fields"][name]
            agg = summ["fields"].setdefault(name, {"n": 0, "sum": 0.0, "l1": 0.0, "l2": 0.0, "linf": 0.0})
            agg["n"] += len(a)
            agg["sum"] += math.fsum(a)
            agg["l1"] += math.fsum(abs(x) for x in a)
            agg["l2"] += math.fsum(x * x for x in a)
            agg["linf"] = max(agg["linf"], max((abs(x) for x in a), default=0.0))
    for agg in summ["fields"].values():
        agg["l2"] = math.sqrt(agg["l2"])
    for g in TEXT_OUTPUT_GLOBS:
        for fn in sorted(glob.glob(os.path.join(rundir, g))):
            t = read_numeric_table(fn)
            if not t:
                continue
            ncol = max(len(r) for r in t)
            cols = [[r[i] for r in t if i < len(r)] for i in range(ncol)]
            summ["text"][os.path.relpath(fn, rundir)] = {
                "rows": len(t),
                "last": t[-1],
                "colsum": [math.fsum(c) for c in cols],
                "colabs": [math.fsum(abs(x) for x in c) for c in cols]}
    return summ


def cmd_bless(args):
    names = [n for n in select_cases(args.cases, args.tags) if os.path.isdir(os.path.join(args.run, n))]
    for n in names:
        rd = os.path.join(args.run, n)
        if run_status(rd).get("status") != "ok":
            print("%-40s skipped (run not ok)" % n)
            continue
        summ = summarize(rd)
        summ["blessed_from"] = run_status(rd).get("reef3d")
        summ["date"] = time.strftime("%Y-%m-%d")
        summ["note"] = args.note or ""
        with open(os.path.join(CASES_DIR, n, "reference.json"), "w") as f:
            json.dump(summ, f, indent=1)
        print("%-40s reference.json written" % n)
    return 0


def rel_close(a, b, rtol, atol):
    return abs(a - b) <= atol + rtol * max(abs(a), abs(b))


def cmd_check(args):
    names = [n for n in select_cases(args.cases, args.tags) if os.path.isdir(os.path.join(args.run, n))]
    bad = 0
    for n in names:
        refp = os.path.join(CASES_DIR, n, "reference.json")
        if not os.path.exists(refp):
            print("%-40s no reference.json" % n)
            continue
        if run_status(os.path.join(args.run, n)).get("status") != "ok":
            print("%-40s FAILED run" % n)
            bad += 1
            continue
        with open(refp) as f:
            ref = json.load(f)
        new = summarize(os.path.join(args.run, n))
        msgs = []
        if ref["final_count"] != new["final_count"]:
            msgs.append("final step %d vs %d" % (ref["final_count"], new["final_count"]))
        for name, r in ref["fields"].items():
            m = new["fields"].get(name)
            if m is None or m["n"] != r["n"]:
                msgs.append("%s: size" % name)
                continue
            scale = r["linf"] * max(r["n"], 1)
            for k in ("l1", "l2", "linf"):
                sc = r["linf"] if k == "linf" else (r[k] if r[k] else scale)
                if abs(m[k] - r[k]) > args.atol + args.rtol * max(abs(sc), abs(m[k])):
                    msgs.append("%s.%s %.6g vs %.6g" % (name, k, m[k], r[k]))
        for fn, r in ref["text"].items():
            m = new["text"].get(fn)
            if m is None or m["rows"] != r["rows"]:
                msgs.append("%s: rows" % fn)
                continue
            for i, (x, y) in enumerate(zip(r["colabs"], m["colabs"])):
                if not rel_close(x, y, args.rtol, args.atol):
                    msgs.append("%s col %d" % (fn, i))
                    break
        print("%-40s %s" % (n, "ok" if not msgs else "DIFF: " + "; ".join(msgs[:6])))
        bad += bool(msgs)
    return 1 if bad else 0


def cmd_list(args):
    for n in select_cases(args.cases, args.tags):
        c = load_case(n)
        print("%-40s np=%d steps=%-4d [%s]\n    %s" % (n, c["np"], c["steps"], ",".join(c["tags"]),
                                                     c.get("description", "")))
    return 0


def cmd_write(args):
    """write the resolved inputs of the cases without running (for inspection / manual runs)."""
    for n in select_cases(args.cases, args.tags):
        write_inputs(load_case(n), os.path.join(args.out, n))
        print(os.path.join(args.out, n))
    return 0


# ----------------------------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def sel(p):
        p.add_argument("--cases", nargs="*", help="case name globs")
        p.add_argument("--tags", nargs="*", help="only cases with one of these tags")

    def runopts(p):
        p.add_argument("--divemesh", required=True)
        p.add_argument("--mpirun", default=os.environ.get("REEF3D_MPIRUN", "mpirun"),
                       help='MPI launcher, e.g. "mpirun --oversubscribe" (env REEF3D_MPIRUN)')
        p.add_argument("--steps", type=int, help="override number of steps for all cases")
        p.add_argument("--timeout", type=int, default=1800)
        p.add_argument("--keep", action="store_true", help="keep grids and bulk output")
        p.add_argument("--max-np", type=int, default=int(os.environ.get("REEF3D_MAX_NP", "0")),
                       help="cap the ranks of every case (A/B on a small machine; env REEF3D_MAX_NP)")

    def tol(p):
        p.add_argument("--rtol", type=float, default=1e-6)
        p.add_argument("--atol", type=float, default=1e-12)
        p.add_argument("--require", choices=LEVELS[:3], default="identical",
                       help="worst acceptable result (default identical = bitwise)")

    p = sub.add_parser("list"); sel(p); p.set_defaults(func=cmd_list)
    p = sub.add_parser("write"); sel(p); p.add_argument("--out", required=True); p.set_defaults(func=cmd_write)
    p = sub.add_parser("run"); sel(p); runopts(p)
    p.add_argument("--reef3d", required=True); p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_run)
    p = sub.add_parser("compare"); p.add_argument("ref"); p.add_argument("new"); tol(p)
    p.add_argument("--cases", nargs="*"); p.add_argument("--report")
    p.set_defaults(func=cmd_compare, cases_list=None)
    p = sub.add_parser("ab"); sel(p); runopts(p); tol(p)
    p.add_argument("--ref-bin", required=True); p.add_argument("--new-bin", required=True)
    p.add_argument("--out", required=True); p.set_defaults(func=cmd_ab)
    p = sub.add_parser("bless"); sel(p); p.add_argument("run"); p.add_argument("--note")
    p.set_defaults(func=cmd_bless)
    p = sub.add_parser("check"); sel(p); p.add_argument("run")
    p.add_argument("--rtol", type=float, default=1e-6); p.add_argument("--atol", type=float, default=1e-10)
    p.set_defaults(func=cmd_check)

    args = ap.parse_args()
    global MAX_NP
    MAX_NP = getattr(args, "max_np", 0) or 0
    sys.exit(args.func(args))


if __name__ == "__main__":
    main()
