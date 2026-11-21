#!/usr/bin/env python3
# Architect: Hans Bihs
"""NTHMP benchmarks 4 and 5: gauges WG1-WG4 against the measurements.

  bm4_eval.py [--bench 5 --h <depth> --name <figure prefix>] <measured file> <out dir> <label>=<run> [...]

Benchmark 4 (default; IRPHE test 17, Grilli et al. 2017; case 28): the gauges 0.6 / 1.6 / 2.6 / 3.6 m from the gate,
the far wall 5.47 m from it, h = 0.33 m.
Benchmark 5 (--bench 5; Viroulet et al. 2014; case 32): the gauges 0.45 / 0.75 / 1.05 / 1.35 m from the gate, the
far wall 1.9 m from it, h from --h; the record starts at the release of the gate (t = 0): the runs are not shifted
(--align crest shifts them as for benchmark 4). "lead": the first crest at WG1 of the run before the measured one.

The measured file is t (s), eta WG1-WG4 (cm). The time of each run is shifted so that its first crest at WG1 (the
first maximum above 20 % of the WG1 maximum) falls on the measured one (the release time of the gate is not in the
record); the same shift for all gauges. Metrics from the first arrival at each gauge to the arrival of the wave
reflected at the far wall of the flume (at sqrt(g h)), at most 4 s.
"""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from slide_eval import style, gauges, slide_log, BLUE, LBLUE, GREEN, ORANGE, GREY

COLS = [BLUE, ORANGE, GREEN, LBLUE]


def first_crest(t, y):
    thr = 0.2*np.max(y)
    k = np.argmax(y > thr)
    while k + 1 < len(y) and y[k + 1] >= y[k]:
        k += 1
    return t[k]


def main():
    a = sys.argv[1:]
    opt = dict(bench="4", h=None, name=None, align=None)
    while a and a[0].startswith("--"):
        opt[a[0][2:]] = a[1]; a = a[2:]
    bm5 = opt["bench"] == "5"
    meas, out = a[0], a[1]
    runs = [x.split("=", 1) for x in a[2:]]
    xg, wall, h = ([0.45, 0.75, 1.05, 1.35], 1.9, 0.148) if bm5 else ([0.6, 1.6, 2.6, 3.6], 5.47, 0.33)
    h = float(opt["h"]) if opt["h"] else h
    name = opt["name"] or ("32_bm5" if bm5 else "28_bm4")
    M = np.loadtxt(meas, comments="#")
    tm = M[:, 0]
    Em = M[:, 1:5] - np.median(M[tm < tm[0] + (0.1 if bm5 else 2.0), 1:5], axis=0)
    # --align crest (default for benchmark 4) / none (default for benchmark 5)
    align = (opt["align"] or ("none" if bm5 else "crest")) == "crest"
    t0m = first_crest(tm, Em[:, 0]) if align else 0.0
    # arrival at each gauge: |eta| above 5 % of its range
    arr = [tm[np.argmax(np.abs(Em[:, k]) > 0.05*np.ptp(Em[:, k]))] for k in range(4)]
    c = np.sqrt(9.81*h)
    tend = [min(arr[k] + 4.0, arr[0] + (2*wall - xg[0] - xg[k])/c) for k in range(4)]
    plt = style()
    fig, ax = plt.subplots(4, 1, figsize=(9, 9), sharex=True)
    res = {}
    for k in range(4):
        ax[k].plot(tm - t0m, Em[:, k], color="k", lw=1.6, label="measured" if bm5 else "measured (test 17)")
    for n, (lab, run) in enumerate(runs):
        xs, G = gauges(run)
        ts = first_crest(G[:, 0], 100*G[:, 1]) if align else 0.0
        # lead: the first crest at WG1 before the measured one (s)
        r = dict(shift=float(t0m - ts), lead=float(first_crest(tm, Em[:, 0]) - first_crest(G[:, 0], 100*G[:, 1])))
        for k in range(4):
            e = 100*G[:, k + 1]
            ax[k].plot(G[:, 0] - ts, e, color=COLS[n % 4], lw=1.2, label=lab)
            w0, w1 = arr[k], tend[k]
            mm = (tm >= w0) & (tm <= w1)
            ei = np.interp(tm[mm] - t0m, G[:, 0] - ts, e)
            em = Em[mm, k]
            r[f"WG{k+1}"] = dict(max=round(float(ei.max()), 3), min=round(float(ei.min()), 3),
                                 max_meas=round(float(em.max()), 3), min_meas=round(float(em.min()), 3),
                                 nrms=round(float(np.sqrt(np.mean((ei - em)**2))/np.ptp(em)), 3))
        L = slide_log(run)
        r["runout_x"] = float(np.max(L[:, 9])); r["slide_umax"] = float(np.max(L[:, 7]))
        res[lab] = r
    for k in range(4):
        ax[k].axvspan(arr[k] - t0m, tend[k] - t0m, color=GREY, alpha=0.08, lw=0)
        ax[k].set_ylabel(f"WG{k+1} $\\eta$ (cm)"); ax[k].axhline(0, color=GREY, lw=0.5)
    ax[0].legend(fontsize=7, ncol=4)
    ax[-1].set_xlabel("t − t(first crest at WG1) (s)" if align else "t (s), from the release of the gate")
    ax[-1].set_xlim(*(((-0.7, 1.9) if align else (0.0, 2.5)) if bm5 else (-1.0, 7.0)))
    fig.tight_layout(); fig.savefig(os.path.join(out, name + "_gauges.png"))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
