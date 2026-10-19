# The Weber junction runs of weber_case.py against the measurements: the water surface along the main channel
# (Huang et al. 2002, Fig. 5, RUN5), its minimum in the separation zone and its rise across the junction (0.95 and
# +5.0 % for RUN5, 0.92 and +8.8 % for RUN2), the separation zone (reversed depth-averaged u) and the surface map.
# The state is the AMR dump of the run (REEF3D_AMR_DUMP=<case>/amr); z/H0 is the depth over the mean depth across
# x* = -4 (the first points of Fig. 5 are 1.000 there).
# usage: python3 weber_compare.py case [case ...] [--np 1] [--png out.png] [--labels a,b]
import argparse, json, os, sys, numpy as np
T = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, T)
from weber_lib import *

ap = argparse.ArgumentParser()
ap.add_argument("cases", nargs="+")
ap.add_argument("--png")
ap.add_argument("--labels")
a = ap.parse_args()

exp = json.load(open(os.path.join(T, "data", "run5_surface_fig5.json")))["profiles"]
R = [load_any(d) for d in a.cases]
labels = a.labels.split(",") if a.labels else [os.path.basename(os.path.normpath(d)) for d in a.cases]
for d, r, lab in zip(a.cases, R, labels):
    main = r["y"] >= 0
    sec = lambda x0: main & (np.abs(r["x"]-x0) < 0.5*r["dx"]+1e-9)
    H0 = r["h"][sec(-4.0 + 0.5*r["dx"]) | sec(-4.0 - 0.5*r["dx"])].mean()
    z = r["h"]/H0
    up = z[sec(1.0 + 0.5*r["dx"]) | sec(1.0 - 0.5*r["dx"])].mean()
    junc = main & (r["x"] < -1.15) & (r["x"] > -4)     # below the corner cell of the branch
    kmin = np.argmin(np.where(junc, z, 9))
    # separation zone: reversed depth-averaged u (u* > 0: towards x* > 0, upstream) in the main channel below x* = -1
    rev = main & (r["x"] < -1) & (r["u"] > 0)
    Ls = (-1 - r["x"][rev].min()) if rev.any() else 0.0
    Ws = r["y"][rev].max() + 0.5*r["dx"] if rev.any() else 0.0
    h7 = r["h"][sec(-7.0 + 0.5*r["dx"]) | sec(-7.0 - 0.5*r["dx"])].mean()
    print("%-10s t %.0f s  H0(x*=-4) %.4f m  h(x*=-7) %.4f  upstream (x*=1) %.4f (+%.1f %%)  min %.4f at x* %.2f y* %.2f  "
          "reversed flow: length %.2f W, width %.2f W" % (lab, r["simtime"], H0, h7, up, 100*(up-1), z[kmin], r["x"][kmin],
                                                       r["y"][kmin], Ls, Ws))
    if r.get("nz", 1) > 1:
        # NHFLOW: the reversed flow in the top and the bottom layer (the separation zone is larger at the surface)
        ext = []
        for lay in (r["nz"]-1, 0):
            rv = main & (r["x"] < -1) & (r["u3"][:, lay] > 0)
            ext.append(((-1 - r["x"][rv].min()) if rv.any() else 0.0, (r["y"][rv].max() + 0.5*r["dx"]) if rv.any() else 0.0))
        print("%-10s reversed flow, top layer: length %.2f W, width %.2f W; bottom layer: length %.2f W, width %.2f W"
              % ("", ext[0][0], ext[0][1], ext[1][0], ext[1][1]))
    r["z"], r["H0"] = z, H0
    if "q* = 0.75" in (r.get("note") or json.load(open(os.path.join(d, "blocks.json")))["note"]):
        err = []
        for k, pts in exp.items():
            yv = float(k)
            for xv, zv in pts:
                err.append(z[at(r, xv, yv)] - zv)
        err = np.array(err)
        print("%-10s Fig. 5 surface: %d points, mean %+.4f, rms %.4f, max |%.4f|" % ("", len(err), err.mean(),
                                                                                  np.sqrt((err**2).mean()), np.abs(err).max()))

if a.png:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    C = ["#2B4C7E", "#4F8FCF", "#E07B39", "#3A9A5B", "#6B6B6B"]
    fig, ax = plt.subplots(len(exp), 1, figsize=(7.5, 10), sharex=True)
    for i, (k, pts) in enumerate(exp.items()):
        yv = float(k)
        p = np.array(pts)
        ax[i].plot(p[:, 0], p[:, 1], "o", mfc="none", color=C[4], ms=4, label="measured (Fig. 5)")
        for m, (r, lab) in enumerate(zip(R, labels)):
            line = (np.abs(r["y"] - yv) < 0.5*r["dx"]+1e-9) & (r["x"] > -4.5) & (r["x"] < 1.5)
            if not line.any():
                yy = r["y"][np.argmin(np.abs(r["y"]-yv))]
                line = (r["y"] == yy) & (r["x"] > -4.5) & (r["x"] < 1.5)
            o = np.argsort(r["x"][line])
            ax[i].plot(r["x"][line][o], r["z"][line][o], "-", color=C[m % 4], lw=1.4, label=lab)
        ax[i].set_ylim(0.92, 1.10)
        ax[i].text(0.98, 0.1, "y/W = %s" % k, transform=ax[i].transAxes, ha="right")
        ax[i].grid(alpha=0.3)
        for s in ("top", "right"):
            ax[i].spines[s].set_visible(False)
        ax[i].axvspan(-1, 0, color="#dddddd", alpha=0.4, lw=0)
    ax[0].legend(loc="upper left", fontsize=8, frameon=False)
    ax[-1].set_xlabel("x/W  (positive upstream, branch at -1 .. 0)")
    ax[len(exp)//2].set_ylabel("z / H0")
    ax[-1].set_xlim(-4.2, 1.2)
    fig.tight_layout()
    fig.savefig(a.png, dpi=150)
    print("wrote", a.png)
