# P8b: the bed of the Nidelva - harbour sediment runs of t2_harbour.py (SED=1): along the river the mean bed change
# of the cross-sections and the difference between the two banks (the transverse bed slope of the bends), the
# bed volume budget of the held inflow reach (S 76), the free river and the harbour blocks, for several runs
# usage: python3 plot_sed.py run [run ...] --labels ... --out prefix [--supply xs xe ys ye]
import sys, os, numpy as np, argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
_argv = sys.argv
sys.argv = sys.argv[:1] + ["x"]
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "plot_harbour.py")).read().split("ap = argparse")[0])
sys.argv = _argv

ap = argparse.ArgumentParser()
ap.add_argument("runs", nargs="+")
ap.add_argument("--labels", nargs="*")
ap.add_argument("--out", default="sed")
ap.add_argument("--y0", type=float, default=-300.0)
ap.add_argument("--supply", nargs=4, type=float, default=[760, 1000, -2200, -1900])
a = ap.parse_args()
labels = a.labels or [os.path.basename(r) for r in a.runs]
cols = [BLUE, ORANGE, GREEN, LBLUE, GREY]

def centres(G):
    X, Y = G["X"], G["Y"]
    xc = 0.25*(X[:-1, :-1]+X[1:, :-1]+X[:-1, 1:]+X[1:, 1:])
    yc = 0.25*(Y[:-1, :-1]+Y[1:, :-1]+Y[:-1, 1:]+Y[1:, 1:])
    ar = 0.5*np.abs((X[1:, 1:]-X[:-1, :-1])*(Y[:-1, 1:]-Y[1:, :-1]) - (X[:-1, 1:]-X[1:, :-1])*(Y[1:, 1:]-Y[:-1, :-1]))
    return xc, yc, ar

def inbox(xc, yc, b):
    return (xc >= b[0]) & (xc <= b[1]) & (yc >= b[2]) & (yc <= b[3])

rows = []
fig, axs = plt.subplots(3, 1, figsize=(9, 8.5), sharex=True)
for k, (d, lab) in enumerate(zip(a.runs, labels)):
    L0a, PBa, nums = snapshot(d)
    L0b, PBb, _ = snapshot(d, n=nums[0])
    s0, ss, mean, bank, spd = None, [], [], [], []
    vol = {"held": 0.0, "river": 0.0, "harbour": 0.0, "harbour_dep": 0.0, "harbour_ero": 0.0, "river_ero": 0.0,
           "river_dep": 0.0}
    P, H = [], []
    for Ga, Gb in zip(L0a, L0b):
        xc, yc, ar = centres(Ga)
        db = Ga["bed"] - Gb["bed"]
        wet = Gb["WL"] > 0.05
        riv = (yc < a.y0)
        held = inbox(xc, yc, a.supply) & riv
        vol["held"] += np.sum((db*ar)[held])
        free = riv & ~held
        vol["river"] += np.sum((db*ar)[free])
        vol["river_ero"] += np.sum((np.minimum(db, 0)*ar)[free])
        vol["river_dep"] += np.sum((np.maximum(db, 0)*ar)[free])
        m = riv.all(1)
        jm = db.shape[1]//2
        P.append(np.c_[xc[m, jm], yc[m, jm]])
        H.append(held[m][:, 1:-1].all(1))
        mean.append(np.array([np.sum((db*ar*wet)[i])/max(np.sum((ar*wet)[i]), 1e-9) for i in np.where(m)[0]]))
        # the difference between the bed changes of the cells next to the two banks (wet cells, j 1 and n-2)
        bank.append(db[m, 1] - db[m, -2])
        spd.append(np.hypot(Ga["u"], Ga["v"])[m, jm])
    for Ga, Gb in zip(PBa, PBb):
        if Ga["bed"].shape != Gb["bed"].shape:
            continue
        xc, yc, ar = centres(Ga)
        db = Ga["bed"] - Gb["bed"]
        vol["harbour"] += np.sum(db*ar)
        vol["harbour_dep"] += np.sum(np.maximum(db, 0)*ar)
        vol["harbour_ero"] += np.sum(np.minimum(db, 0)*ar)
    P = np.concatenate(P)
    s = np.r_[0, np.cumsum(np.hypot(*np.diff(P, axis=0).T))]
    H = np.concatenate(H)
    shold = s[H].max() if H.any() else 0.0
    axs[0].plot(s, np.concatenate(spd), color=cols[k % 5], lw=1, label=lab)
    axs[1].plot(s, np.concatenate(mean), color=cols[k % 5], lw=1, label=lab)
    axs[2].plot(s, np.concatenate(bank), color=cols[k % 5], lw=1, label=lab)
    rows.append((lab, L0a[0]["time"], L0b[0]["time"], vol))
for ax in axs:
    ax.axvspan(0, shold, color=GREY, alpha=0.12, lw=0)
axs[0].set_ylabel("speed, middle (m/s)")
axs[1].set_ylabel("mean bed change (m)")
axs[2].set_ylabel("bed change, right - left bank (m)")
axs[2].set_xlabel("distance along the river from the inflow (m); grey: the cross-sections held in full (S 76 box of --supply)")
axs[0].legend(fontsize=8)
axs[0].set_title("Nidelva: the bed along the river corridor (block 0)")
fig.tight_layout()
fig.savefig("%s_profile.png" % a.out, dpi=150)
plt.close(fig)

with open("%s_budget.txt" % a.out, "w") as f:
    f.write("# bed volume change (m3) between the first and the last output\n")
    f.write("# run, t0, t1, held reach, free river (net, erosion, deposition), harbour blocks (net, erosion, deposition)\n")
    for lab, t1, t0, v in rows:
        line = "%-22s %6.0f %6.0f  held %9.3g  river %10.4g %10.4g %10.4g  harbour %10.4g %10.4g %10.4g" % (
            lab, t0, t1, v["held"], v["river"], v["river_ero"], v["river_dep"], v["harbour"], v["harbour_ero"],
            v["harbour_dep"])
        f.write(line + "\n")
        print(line)
print("written", a.out)
