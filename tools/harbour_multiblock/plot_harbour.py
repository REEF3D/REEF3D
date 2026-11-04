# P5d: figures of the Nidelva - harbour runs of t2_harbour.py (SFLOW: REEF3D_SFLOW_AMR/*.vts of level 0 (the river
# corridor, block 0) and of the blocks 1.. and patches; NHFLOW: REEF3D_NHFLOW_AMR/*.vts as written by nhflow_amr)
# usage: python3 plot_harbour.py rundir [rundir2 ...] --out prefix
import sys, os, re, glob, argparse, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

BLUE, LBLUE, ORANGE, GREEN, GREY = "#2B4C7E", "#4F8FCF", "#E07B39", "#3A9A5B", "#6B6B6B"
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
                     "grid.alpha": 0.3})

WD = 10.0

def readvts(fn):
    s = open(fn).read()
    ext = [int(v) for v in re.search(r'WholeExtent="([^"]+)"', s).group(1).split()]
    nx, ny = ext[1], ext[3]
    out = {}
    for m in re.finditer(r'<DataArray type="Float(?:64|32)" Name="(\w+)"[^>]*format="ascii">(.*?)</DataArray>', s, re.S):
        a = np.array(m.group(2).split(), float)
        if a.size == nx*ny:
            out[m.group(1)] = a.reshape(ny, nx).T
    m = re.search(r'<Points>\s*<DataArray[^>]*>(.*?)</DataArray>', s, re.S)
    P = np.array(m.group(1).split(), float).reshape(ny+1, nx+1, 3).transpose(1, 0, 2)
    t = re.search(r'Name="TimeValue"[^>]*>\s*([-0-9.eE+]+)', s)
    lv = re.search(r'Name="level"[^>]*>\s*(-?\d+)', s)
    out["X"], out["Y"] = P[:, :, 0], P[:, :, 1]
    if "u" not in out and "u_top" in out:          # NHFLOW: the velocity of the top layer, WL from eta and the bed
        out["u"], out["v"] = out["u_top"], out["v_top"]
    if "WL" not in out:
        out["WL"] = out["eta"] + WD - out["bed"]
    out["time"] = float(t.group(1)) if t else np.nan
    out["level"] = int(lv.group(1)) if lv else 0
    return out

def out_time(f):
    m = re.search(r'Name="TimeValue"[^>]*>\s*([-0-9.eE+]+)', open(f).read(2000))
    return float(m.group(1)) if m else np.nan

def snapshot(d, n=None, t=None):
    """the grids of output n (default: the last, or the one nearest to the time t): level 0 (the corridor) and the
    blocks 1.. / patches"""
    for sub in ("REEF3D_SFLOW_AMR", "REEF3D_NHFLOW_AMR"):
        if os.path.isdir(os.path.join(d, sub)):
            break
    files = sorted(glob.glob(os.path.join(d, sub, "*-L0-*.vts")))
    nums = sorted({int(re.search(r"-L0-(\d+)-", f).group(1)) for f in files})
    if t is not None:
        n = min(nums, key=lambda q: abs(out_time(sorted(glob.glob(os.path.join(d, sub, "*-L0-%08d-*.vts" % q)))[0])-t))
    n = nums[-1] if n is None else n
    L0 = [readvts(f) for f in sorted(glob.glob(os.path.join(d, sub, "*-L0-%08d-*.vts" % n)))]
    PB = [readvts(f) for f in sorted(glob.glob(os.path.join(d, sub, "*-%08d-*-*.vts" % n))) if "-L0-" not in f]
    return L0, PB, nums

def cells(G, f, y0=None, wet=0.05):
    """polygons and values of the wet cells (level-0 cells of the corridor in the harbour block, the ring, left out)"""
    X, Y = G["X"], G["Y"]
    poly = np.stack([np.stack([X[:-1, :-1], Y[:-1, :-1]], -1), np.stack([X[1:, :-1], Y[1:, :-1]], -1),
                     np.stack([X[1:, 1:], Y[1:, 1:]], -1), np.stack([X[:-1, 1:], Y[:-1, 1:]], -1)], 2)
    m = G["WL"] > wet
    if y0 is not None:
        m &= 0.25*(Y[:-1, :-1]+Y[1:, :-1]+Y[:-1, 1:]+Y[1:, 1:]) < y0
    v = f(G)
    return poly[m], v[m]

def speed(G):
    return np.hypot(G["u"], G["v"])

def mapplot(ax, grids, f, cmap, vmin, vmax, y0):
    pc = None
    for G, cut in grids:
        p, v = cells(G, f, y0 if cut else None)
        pc = PolyCollection(p, array=v, cmap=cmap, edgecolors="none")
        pc.set_clim(vmin, vmax)
        ax.add_collection(pc)
    ax.set_aspect("equal")
    return pc

ap = argparse.ArgumentParser()
ap.add_argument("runs", nargs="+")
ap.add_argument("--labels", nargs="*")
ap.add_argument("--out", default="harbour")
ap.add_argument("--y0", type=float, default=-300.0)
ap.add_argument("--wd", type=float, default=10.0, help="still water level (NHFLOW output has no WL)")
ap.add_argument("--vmax", type=float, default=1.5)
ap.add_argument("--time", type=float, help="the output nearest to this time (default: the last)")
a = ap.parse_args()
WD = a.wd
labels = a.labels or [os.path.basename(r) for r in a.runs]

# ---- speed and water level of the last output of every run
for d, lab in zip(a.runs, labels):
    L0, PB, nums = snapshot(d, t=a.time)
    grids = [(G, True) for G in L0] + [(G, False) for G in PB]
    fig, axs = plt.subplots(1, 2, figsize=(11, 6.5))
    pc = mapplot(axs[0], grids, speed, "viridis", 0.0, a.vmax, a.y0)
    fig.colorbar(pc, ax=axs[0], shrink=0.7, label="speed (m/s)")
    eta0 = np.nanmedian(np.concatenate([G["eta"][G["WL"] > 0.05] for G in PB]))
    pc = mapplot(axs[1], grids, lambda G: G["eta"], "RdBu_r", None, None, a.y0)
    e = np.concatenate([cells(G, lambda g: g["eta"], a.y0 if cut else None)[1] for G, cut in grids])
    pc.set_clim(np.percentile(e, 1), np.percentile(e, 99))
    for c in axs[1].collections:
        c.set_clim(np.percentile(e, 1), np.percentile(e, 99))
    fig.colorbar(pc, ax=axs[1], shrink=0.7, label="eta (m)")
    for ax, t in zip(axs, ("speed", "water level")):
        ax.set_xlim(-450, 2350); ax.set_ylim(-2250, 2050)
        ax.set_title("%s, %s at t = %.0f s" % (lab, t, L0[0]["time"]))
        ax.set_xlabel("x (m)"); ax.set_ylabel("y (m)")
    fig.tight_layout()
    fig.savefig("%s_%s_map.png" % (a.out, lab), dpi=150)
    plt.close(fig)

    # zoom: the mouth, the canal and the piers
    fig, ax = plt.subplots(figsize=(8, 6))
    pc = mapplot(ax, grids, speed, "viridis", 0.0, a.vmax, a.y0)
    for G, cut in grids:
        if not cut:
            ax.quiver(*[0.25*(G[k][:-1, :-1]+G[k][1:, :-1]+G[k][:-1, 1:]+G[k][1:, 1:])[::4, ::4] for k in ("X", "Y")],
                      np.where(G["WL"] > 0.05, G["u"], np.nan)[::4, ::4], np.where(G["WL"] > 0.05, G["v"], np.nan)[::4, ::4],
                      scale=25, width=0.0015, color="w")
    fig.colorbar(pc, ax=ax, shrink=0.8, label="speed (m/s)")
    ax.set_xlim(-300, 2000); ax.set_ylim(-500, 1500)
    ax.set_title("%s: the river mouth, the canal and the harbour basins (t = %.0f s)" % (lab, L0[0]["time"]))
    fig.tight_layout()
    fig.savefig("%s_%s_zoom.png" % (a.out, lab), dpi=150)
    plt.close(fig)

# ---- along the river: water level and speed on the middle line of the corridor (level 0), all runs
fig, axs = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
cols = [BLUE, ORANGE, GREEN, LBLUE, GREY]
for k, (d, lab) in enumerate(zip(a.runs, labels)):
    L0, PB, nums = snapshot(d, t=a.time)
    xs, et, sp = [], [], []
    for G in L0:
        jm = G["eta"].shape[1]//2
        xc = 0.25*(G["X"][:-1, jm]+G["X"][1:, jm]+G["X"][:-1, jm+1]+G["X"][1:, jm+1])
        yc = 0.25*(G["Y"][:-1, jm]+G["Y"][1:, jm]+G["Y"][:-1, jm+1]+G["Y"][1:, jm+1])
        m = yc < a.y0
        xs.append(np.c_[xc[m], yc[m]]); et.append(G["eta"][m, jm]); sp.append(np.hypot(G["u"], G["v"])[m, jm])
    P = np.concatenate(xs)
    s = np.r_[0, np.cumsum(np.hypot(*np.diff(P, axis=0).T))]
    axs[0].plot(s, np.concatenate(et), color=cols[k % 5], label=lab)
    axs[1].plot(s, np.concatenate(sp), color=cols[k % 5], label=lab)
axs[0].set_ylabel("eta (m)"); axs[1].set_ylabel("speed (m/s)"); axs[1].set_xlabel("distance along the river from the inflow (m)")
axs[0].legend()
axs[0].set_title("Nidelva: the middle line of the river corridor (block 0)")
fig.tight_layout()
fig.savefig("%s_profile.png" % a.out, dpi=150)
plt.close(fig)

# ---- water volume and inflow over time (logs)
fig, axs = plt.subplots(1, 2, figsize=(10, 3.5))
for k, (d, lab) in enumerate(zip(a.runs, labels)):
    lg = glob.glob(os.path.join(d, "REEF3D_*_AMR", "REEF3D_*_AMR_log.dat"))
    if lg:
        A = np.loadtxt(lg[0], comments="#", usecols=(1, 5), ndmin=2)
        axs[0].plot(A[:, 0], (A[:, 1]/A[0, 1]-1.0)*100.0, color=cols[k % 5], label=lab)
    t, q = [], []
    for l in open(os.path.join(d, "reef.log")):
        if l.startswith("simtime:"):
            tt = float(l.split()[1])
        if l.startswith("Inflow:"):
            t.append(tt); q.append(float(l.split()[1]))
    axs[1].plot(t, q, color=cols[k % 5], label=lab)
axs[0].set_xlabel("t (s)"); axs[0].set_ylabel("water volume change (%)")
axs[0].set_title("water volume (all blocks)")
axs[1].set_title("discharge of the native inflow (W 10)")
axs[1].set_xlabel("t (s)"); axs[1].set_ylabel("discharge at the inflow (m3/s)")
axs[0].legend()
fig.tight_layout()
fig.savefig("%s_logs.png" % a.out, dpi=150)
plt.close(fig)
print("written", a.out)
