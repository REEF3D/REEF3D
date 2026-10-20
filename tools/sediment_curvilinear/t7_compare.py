# P7: the runs of t7_talmon_bend.py (runs/t7_talmon<TAG>/t7.npz) against the measurements of Talmon & de Graaff (1989)
# run 3 and against each other: per run the rms of a/a0 (all, sections 6-24, 25-45), outer - inner, the deepest outer
# point; a figure of the depth at the outer bank (0.45 m from the inner bank), outer - inner and the inner bank along
# the 48 sections, one column per grid.
# usage: python3 t7_compare.py --grids _ds05,_ds025 --runs ",_s952,_lag" --labels "S 95 1,S 95 2,S 95 2 + S 64 0.75" [--png out.png]
import os, json, argparse, numpy as np
T = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser()
ap.add_argument("--grids", default="_ds05,_ds025")
ap.add_argument("--glabels", default="10 cells across (ds 0.05 m),20 cells across (ds 0.025 m)")
ap.add_argument("--runs", default=",_s952,_lag")
ap.add_argument("--labels", default="S 95 1 (P6c),S 95 2 (spiral lag),S 95 2 + S 64 0.75 (spiral and suspension lags)")
ap.add_argument("--png")
a = ap.parse_args()
grids, glabels = a.grids.split(","), a.glabels.split(",")
runs, labels = a.runs.split(","), a.labels.split(",")
D = json.load(open(os.path.join(T, "data", "talmon1989_run3.json")))
A_mea = np.array(D["a_over_a0"])
ks_ = np.arange(1, 49)
rms = lambda q: np.sqrt(np.mean(q**2))
res = {}
for g in grids:
    for r in runs:
        f = os.path.join(T, "runs", "t7_talmon%s%s" % (g, r), "t7.npz")
        if not os.path.exists(f):
            continue
        A = np.load(f)["A_mod"]
        res[(g, r)] = A
        dd = A - A_mea
        do_m, do_s = A_mea[:, -1] - A_mea[:, 0], A[:, -1] - A[:, 0]
        kp = np.argmax(A[:, -1])
        print("%-7s %-6s rms all %.3f  6-24 %.3f  25-45 %.3f | outer-inner 12-18 %.2f (meas. %.2f)  25-45 %.2f (%.2f) | deepest outer section %d %.2f a0 (meas. 16 1.73)"
              % (g, r or "_s951", rms(dd), rms(dd[5:24]), rms(dd[24:45]), do_s[11:18].mean(), do_m[11:18].mean(),
                 do_s[24:45].mean(), do_m[24:45].mean(), kp+1, A[kp, -1]))
for r in runs:
    if all((g, r) in res for g in grids[:2]) and len(grids) > 1:
        A1, A2 = res[(grids[0], r)], res[(grids[1], r)]
        print("grids %s vs %s, %-6s: rms of the difference %.3f a0, outer bank %.3f a0" % (grids[0], grids[1], r or "_s951", rms(A1 - A2), rms(A1[:, -1] - A2[:, -1])))
if a.png:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True, "grid.color": "#e6e6e6"})
    cols = ["#6B6B6B", "#4F8FCF", "#E07B39", "#3A9A5B", "#2B4C7E"]
    if len(runs) <= 3:
        cols = ["#6B6B6B", "#4F8FCF", "#2B4C7E"]
    fig, axs = plt.subplots(3, len(grids), figsize=(max(5.6*len(grids), 8.0), 9.0), sharex=True, sharey="row", squeeze=False)
    rows = [("a/a0 at the outer bank (0.45 m)", lambda A: A[:, -1]),
            ("(a outer - a inner)/a0", lambda A: A[:, -1] - A[:, 0]),
            ("a/a0 at the inner bank (0.05 m)", lambda A: A[:, 0])]
    for ic, (g, gl) in enumerate(zip(grids, glabels)):
        for ir, (yl, fn) in enumerate(rows):
            ax = axs[ir, ic]
            ax.plot(ks_, fn(A_mea), "o", ms=3.2, mfc="none", mec="#222222", mew=0.9, label="measured")
            for k, (r, lab) in enumerate(zip(runs, labels)):
                if (g, r) in res:
                    last = (k == len(runs) - 1)
                    ax.plot(ks_, fn(res[(g, r)]), color=cols[k % len(cols)], lw=2.2 if last else (1.1 if k == 0 else 1.3),
                            ls="--" if k == 0 else "-", label=lab, zorder=3 if last else 2)
            ax.axvline(4.75, color="#6B6B6B", lw=0.7, ls=":"); ax.axvline(45, color="#6B6B6B", lw=0.7, ls=":")
            if ic == 0:
                ax.set_ylabel(yl)
            if ir == 0:
                ax.set_title(gl, fontsize=9)
            if ir == 2:
                ax.set_xlabel("cross-section (bend from 4.75 to 45)")
    h_, l_ = axs[0, 0].get_legend_handles_labels()
    fig.legend(h_, l_, loc="lower center", ncol=2, frameon=False, fontsize=8)
    fig.suptitle("Talmon & de Graaff (1989) run 3, LFM curved flume: SFLOW with the lags of P7", fontsize=10)
    fig.tight_layout(rect=(0, 0.035*((len(l_)+1)//2), 1, 0.97))
    fig.savefig(a.png, dpi=150)
    print("wrote", a.png)
