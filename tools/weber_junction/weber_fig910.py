# Vertical velocity profiles of NHFLOW runs of weber_case.py against the measurements of Huang et al. (2002), Figs. 9
# (streamwise: -u/U0 in the main channel, v/U0 in the branch) and 10 (transverse: v/U0 in the main channel, -u/U0 in
# the branch), RUN5 (data/fig9_10_profiles.json, digitized by digitize_fig9_10.py).  The model per sigma layer at the
# station (linear in x and y between the cell centres), then linear in z.  Printed per station: the rms and the mean of
# model - measured; per figure the rms over all points.
# usage: python3 weber_fig910.py case [case ...] [--labels a,b] [--png out.png]
import argparse, json, os, sys, numpy as np
from scipy.interpolate import LinearNDInterpolator
T = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, T)
from weber_lib import *

def profiles(r, stations):
    """model profile (z/W of the layers, value/U0) at every station"""
    pts = np.c_[r["x"], r["y"]]
    vals = np.concatenate([r["u3"], r["v3"], r["h"][:, None]], axis=1)
    near = lambda s: (np.abs(r["x"] - s["x"]) < 1.5) & (np.abs(r["y"] - s["y"]) < 1.5)
    out = []
    for s in stations:
        k = near(s)
        f = LinearNDInterpolator(pts[k], vals[k])
        q = f(np.array([[s["x"], s["y"]]]))[0]
        nz = r["nz"]
        u3, v3, h = q[:nz], q[nz:2*nz], q[-1]
        zc = (np.arange(nz) + 0.5)/nz*h/W
        val = {"-u": -u3, "v": v3}[s["q"]]/U0
        out.append((zc, val, h/W))
    return out

def score(stations, prof):
    res = []
    for s, (zc, val, hw) in zip(stations, prof):
        z = np.array(s["z"]); m = np.array(s["val"])
        mod = np.interp(z, zc, val)
        res.append(dict(d=mod - m, n=len(z)))
    return res

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cases", nargs="+")
    ap.add_argument("--labels")
    ap.add_argument("--png")
    a = ap.parse_args()
    data = json.load(open(os.path.join(T, "data", "fig9_10_profiles.json")))
    st = data["stations"]
    labels = a.labels.split(",") if a.labels else [os.path.basename(c.rstrip("/")) for c in a.cases]
    runs = []
    for c, lab in zip(a.cases, labels):
        r = load_any(c)
        pr = profiles(r, st)
        sc = score(st, pr)
        runs.append((lab, pr, sc))
        print("%-12s t %.0f s" % (lab, r["simtime"]))
        for fig in ("fig9", "fig10"):
            dd = np.concatenate([x["d"] for s, x in zip(st, sc) if s["fig"] == fig])
            print("   %-5s %3d points  rms %.3f  mean %+.3f" % (fig, len(dd), np.sqrt(np.mean(dd**2)), dd.mean()))
            line = []
            for s, x in zip(st, sc):
                if s["fig"] != fig:
                    continue
                line.append("(%s %+.2f,%.2f) %+.2f/%.2f" % (s["panel"], s["x"], s["y"], x["d"].mean(), np.sqrt(np.mean(x["d"]**2))))
            for i in range(0, len(line), 3):
                print("        " + "   ".join(line[i:i+3]))
    if a.png:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.spines.top": False,
                             "axes.spines.right": False, "axes.grid": True, "grid.color": "#e6e6e6"})
        cols = ["#2B4C7E", "#E07B39", "#3A9A5B", "#4F8FCF"]
        # a station per panel: Fig. 9 at x/W 1, -1, -2, -6 and Fig. 10 at x/W -1, -2 (y/W 0.25, 0.5, 0.75)
        rows = [("fig9", (1.0, -1.0)), ("fig9", (-2.0, -6.0)), ("fig10", (-1.0, -2.0))]
        fig, axs = plt.subplots(3, 6, figsize=(14, 9.5), sharey=True)
        for ir, (fg, xs_) in enumerate(rows):
            for ic in range(6):
                ax = axs[ir, ic]
                xq, yq = xs_[ic//3], (0.25, 0.5, 0.75)[ic % 3]
                ii = [i for i, s in enumerate(st) if s["fig"] == fg and s["x"] == xq and s["y"] == yq]
                if not ii:
                    ax.set_visible(False); continue
                i = ii[0]; s = st[i]
                ax.plot(s["val"], s["z"], "o", ms=3.2, mfc="none", mec="#6B6B6B", mew=0.9, label="measured")
                for k, (lab, pr, _) in enumerate(runs):
                    zc, val, hw = pr[i]
                    ax.plot(np.r_[val, val[-1]], np.r_[zc, hw], color=cols[k % 4], lw=1.4, label=lab)
                ax.axvline(0, color="#6B6B6B", lw=0.6)
                ax.set_title("Fig. %s  x/W = %g, y/W = %g" % (fg[3:], xq, yq), fontsize=8)
                ax.set_xlabel(("-u/U0 (streamwise)" if s["q"] == "-u" else "v/U0 (transverse)"), fontsize=7.5)
                ax.set_ylim(0, 0.37)
                if fg == "fig9":
                    ax.set_xlim(0.3, 1.6)
                else:
                    ax.set_xlim(-0.35, 0.6)
            axs[ir, 0].set_ylabel("z/W")
        h_, l_ = axs[0, 0].get_legend_handles_labels()
        fig.legend(h_, l_, loc="lower center", ncol=len(l_), frameon=False, fontsize=8)
        fig.suptitle("RUN5 velocity profiles: measured (Shumate 1998, in Huang et al. 2002, Figs. 9 and 10) and NHFLOW", fontsize=10)
        fig.tight_layout(rect=(0, 0.03, 1, 1))
        fig.savefig(a.png, dpi=150)
        print("wrote", a.png)
