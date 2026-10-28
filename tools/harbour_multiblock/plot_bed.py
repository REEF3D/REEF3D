# P5d: bed change of the SFLOW sediment run of t2_harbour.py (SED=1): bed at the last output minus the first, the
# river corridor (block 0) and the harbour blocks
# usage: python3 plot_bed.py rundir out.png [vmax [xs xe ys ye]]   (the window of the right panel, default the loop)
import sys, os, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
sys.argv, args = sys.argv[:1] + ["x"], sys.argv[1:]
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "plot_harbour.py")).read().split("ap = argparse")[0])
d, out = args[0], args[1]
vmax = float(args[2]) if len(args) > 2 else 0.5
zoom = [float(v) for v in args[3:7]] if len(args) > 6 else [-350, 1000, -1050, -300]

L0a, PBa, nums = snapshot(d, n=None)
L0b, PBb, _ = snapshot(d, n=nums[0])
fig, axs = plt.subplots(1, 2, figsize=(13, 6.5), gridspec_kw={"width_ratios": [1, 1.4]})
for ax, lim in zip(axs, ([-450, 2350, -2250, 2050], zoom)):
    pc = None
    for Ga, Gb, cut in [(a_, b_, True) for a_, b_ in zip(L0a, L0b)] + [(a_, b_, False) for a_, b_ in zip(PBa, PBb)]:
        if Ga["bed"].shape != Gb["bed"].shape:
            continue
        G = dict(Ga)
        G["db"] = Ga["bed"] - Gb["bed"]
        p, v = cells(G, lambda g: g["db"], -300.0 if cut else None, wet=0.05)
        pc = PolyCollection(p, array=v, cmap="RdBu_r", edgecolors="none")
        pc.set_clim(-vmax, vmax)
        ax.add_collection(pc)
    ax.set_xlim(lim[0], lim[1]); ax.set_ylim(lim[2], lim[3]); ax.set_aspect("equal")
    ax.set_xlabel("x (m)"); ax.set_ylabel("y (m)")
fig.colorbar(pc, ax=axs, shrink=0.7, label="bed change (m)")
axs[0].set_title("bed change after %.0f s (S 35 morph. factor)" % L0a[0]["time"])
axs[1].set_title("x %.0f..%.0f, y %.0f..%.0f m" % tuple(zoom) if len(args) > 6 else "the loop of the corridor (S 86 bed slope, S 95 spiral flow)")
fig.savefig(out, dpi=150, bbox_inches="tight")
print("written", out)
