#!/usr/bin/env python3
# REEF3D::CORAL against MoorDyn v2: figures of the time series of coral_moordyn_compare (out/<case>_coral.csv,
# out/<case>_moordyn.csv), from Regression/moordyn:
#   ./plot_compare.py [--out <folder>] [case ...]
# One figure per case (out/fig/<case>.png: the force on the body or its motion, the end tensions of the lines)
# and the overview out/fig/overview.png (one channel per case). CORAL solid blue, MoorDyn dashed orange.
# Architect: Hans Bihs
import argparse
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLUE, LBLUE, ORANGE, GREEN, GREY = "#2B4C7E", "#4F8FCF", "#E07B39", "#3A9A5B", "#6B6B6B"
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Liberation Sans", "Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.color": "#DDDDDD",
    "grid.linewidth": 0.6,
    "axes.edgecolor": GREY,
    "axes.labelcolor": "#333333",
    "xtick.color": GREY,
    "ytick.color": GREY,
    "legend.frameon": False,
})

# the channel of each case for the overview, its label and scale
KEY = {
    "oc3_static": ("L1_Tb", "fairlead tension [kN]", 1e-3),
    "oc3_surge": ("L1_Tb", "fairlead tension [kN]", 1e-3),
    "oc3_heave": ("L1_Tb", "fairlead tension [kN]", 1e-3),
    "oc3_zeta": ("L1_Tb", "fairlead tension [kN]", 1e-3),
    "taut_poly": ("L1_Tb", "fairlead tension [kN]", 1e-3),
    "poly_curve": ("L1_Tb", "fairlead tension [kN]", 1e-3),
    "chain_snap": ("L1_Tb", "fairlead tension [kN]", 1e-3),
    "chain_surface": ("L1_Tb", "fairlead tension [kN]", 1e-3),
    "oc3_current": ("Fy", "cross force on the body [kN]", 1e-3),
    "oc3_shear": ("L1_Tb", "fairlead tension [kN]", 1e-3),
    "bridle": ("M1_Ta", "chain tension at the bridle point [kN]", 1e-3),
    "clump": ("F1_Tb", "fairlead tension [kN]", 1e-3),
    "buoy": ("U1_z", "buoy height [m]", 1.0),
    "chain_wire_chain": ("S3_Tb", "fairlead tension [kN]", 1e-3),
    "spar_decay": ("x", "surge [m]", 1.0),
    "spar_forced": ("x", "surge [m]", 1.0),
}


def flat_ylim(ax, *series):
    """a channel that hardly varies (a static case): +-0.5 % about its mean instead of the round-off"""
    v = [x for s in series for x in s]
    lo, hi = min(v), max(v)
    mid = 0.5*(lo + hi)
    if hi - lo < 1e-4*max(abs(mid), 1e-12):
        d = max(5e-3*abs(mid), 1e-6)
        ax.set_ylim(mid - d, mid + d)
        ax.ticklabel_format(axis="y", useOffset=False)


def read(path):
    with open(path) as f:
        rows = list(csv.reader(f))
    names = rows[0]
    data = {n: [] for n in names}
    for r in rows[1:]:
        for n, v in zip(names, r):
            data[n].append(float(v))
    return data


def channels(d):
    """the channels of a case figure: motion of a free body or the force on the body, then the end tensions"""
    ch = []
    if max(abs(v) for v in d["x"]) > 0 and "spar" in d.get("_name", ""):
        ch += [("x", "surge [m]", 1.0), ("pitch", "pitch [deg]", 57.29578)]
    ch += [("Fx", "Fx [kN]", 1e-3), ("Fz", "Fz [kN]", 1e-3)]
    for n in d:
        if n.endswith("_Ta") or n.endswith("_Tb"):
            ch.append((n, n.replace("_", " ") + " [kN]", 1e-3))
    return ch[:8]


def case_figure(name, c, m, outdir):
    c["_name"] = name
    ch = channels(c)
    n = len(ch)
    ncol = 2
    nrow = (n + 1)//2
    fig, axs = plt.subplots(nrow, ncol, figsize=(10, 2.2*nrow), sharex=True, squeeze=False)
    for k, (key, label, s) in enumerate(ch):
        ax = axs[k//ncol][k % ncol]
        ax.plot(c["t"], [v*s for v in c[key]], color=BLUE, lw=1.4, label="CORAL")
        ax.plot(m["t"], [v*s for v in m[key]], color=ORANGE, lw=1.2, ls="--", label="MoorDyn 2.7.2")
        flat_ylim(ax, [v*s for v in c[key]], [v*s for v in m[key]])
        ax.set_ylabel(label)
        if k == 0:
            ax.legend(loc="best")
    for k in range(n, nrow*ncol):
        axs[k//ncol][k % ncol].axis("off")
    for ax in axs[-1]:
        ax.set_xlabel("t [s]")
    fig.suptitle(name.replace("_", " "), x=0.02, ha="left", color="#333333")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, name + ".png"), dpi=150)
    plt.close(fig)


def overview(cases, outdir, src):
    n = len(cases)
    ncol = 4
    nrow = (n + ncol - 1)//ncol
    fig, axs = plt.subplots(nrow, ncol, figsize=(16, 2.6*nrow), squeeze=False)
    for k, name in enumerate(cases):
        ax = axs[k//ncol][k % ncol]
        c = read(os.path.join(src, name + "_coral.csv"))
        m = read(os.path.join(src, name + "_moordyn.csv"))
        key, label, s = KEY.get(name, ("Fx", "Fx [kN]", 1e-3))
        ax.plot(c["t"], [v*s for v in c[key]], color=BLUE, lw=1.3, label="CORAL")
        ax.plot(m["t"], [v*s for v in m[key]], color=ORANGE, lw=1.1, ls="--", label="MoorDyn")
        flat_ylim(ax, [v*s for v in c[key]], [v*s for v in m[key]])
        ax.set_title(name.replace("_", " "), fontsize=10, loc="left", color="#333333")
        ax.set_ylabel(label, fontsize=8)
        ax.tick_params(labelsize=8)
        if k == 0:
            ax.legend(loc="best", fontsize=8)
    for k in range(n, nrow*ncol):
        axs[k//ncol][k % ncol].axis("off")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "overview.png"), dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="out")
    ap.add_argument("cases", nargs="*")
    a = ap.parse_args()
    names = a.cases or [f[:-len("_coral.csv")] for f in sorted(os.listdir(a.out)) if f.endswith("_coral.csv")]
    names = [n for n in names if os.path.exists(os.path.join(a.out, n + "_moordyn.csv"))]
    if not names:
        sys.exit("no case results in " + a.out)
    order = list(KEY)
    names.sort(key=lambda n: order.index(n) if n in order else len(order))
    fig_dir = os.path.join(a.out, "fig")
    os.makedirs(fig_dir, exist_ok=True)
    for n in names:
        case_figure(n, read(os.path.join(a.out, n + "_coral.csv")), read(os.path.join(a.out, n + "_moordyn.csv")), fig_dir)
    overview(names, fig_dir, a.out)
    print("figures in " + fig_dir)


if __name__ == "__main__":
    main()
