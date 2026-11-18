#!/usr/bin/env python3
# Architect: Hans Bihs
"""NTHMP benchmark 4 (IRPHE test 17, Grilli et al. 2017): gauges WG1-WG4 of case 28 against the measurements.

  bm4_eval.py <measured file> <out dir> <label>=<run> [<label>=<run> ...]

The measured file is t (s), eta WG1-WG4 (cm). The time of each run is shifted so that its first crest at WG1 (the
first maximum above 20 % of the WG1 maximum) falls on the measured one (the release time of the gate is not in the
record); the same shift for all gauges. Metrics from the first arrival at each gauge to the arrival of the wave
reflected at the far wall of the flume (x = 5.47 m from the gate, at sqrt(g h)), at most 4 s.
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
    meas, out = sys.argv[1], sys.argv[2]
    runs = [a.split("=", 1) for a in sys.argv[3:]]
    M = np.loadtxt(meas, comments="#")
    tm = M[:, 0]; Em = M[:, 1:5] - np.median(M[tm < tm[0] + 2.0, 1:5], axis=0)
    t0m = first_crest(tm, Em[:, 0])
    # arrival at each gauge: |eta| above 5 % of its range
    arr = [tm[np.argmax(np.abs(Em[:, k]) > 0.05*np.ptp(Em[:, k]))] for k in range(4)]
    xg = [0.6, 1.6, 2.6, 3.6]; c = np.sqrt(9.81*0.33)
    tend = [min(arr[k] + 4.0, arr[0] + (2*5.47 - xg[0] - xg[k])/c) for k in range(4)]
    plt = style()
    fig, ax = plt.subplots(4, 1, figsize=(9, 9), sharex=True)
    res = {}
    for k in range(4):
        ax[k].plot(tm - t0m, Em[:, k], color="k", lw=1.6, label="measured (test 17)")
    for n, (lab, run) in enumerate(runs):
        xs, G = gauges(run)
        ts = first_crest(G[:, 0], 100*G[:, 1])
        r = dict(shift=float(t0m - ts))
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
    ax[0].legend(fontsize=7, ncol=4); ax[-1].set_xlabel("t − t(first crest at WG1) (s)")
    ax[-1].set_xlim(-1.0, 7.0)
    fig.tight_layout(); fig.savefig(os.path.join(out, "28_bm4_gauges.png"))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
