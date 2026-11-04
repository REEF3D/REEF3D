# The streamwise velocity in the cross section x/W = -2 (one W below the junction) of NHFLOW runs of weber_case.py
# against the measurements of Huang et al. (2002), Fig. 8 (data/fig8_u_xm2.json, digitized from the colours: bins of 0.1,
# saturated at -1.4).  u/U0 positive upstream, U0 = 0.628 m/s.  The model per sigma layer, interpolated to the measured
# points (linear in y and in z); the model clipped to the colour bar (-1.4 .. 0.1) for the rms.  The width of the
# separation zone: the first y/W from the branch-side wall with u/U0 <= -0.5, near the bed and near the surface.
# usage: python3 weber_fig8.py --run5 case [--run2 case] [--labels a] [--png out.png]
import argparse, json, os, sys, numpy as np
T = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, T)
from weber_lib import *

ap = argparse.ArgumentParser()
ap.add_argument("--run5")
ap.add_argument("--run2")
ap.add_argument("--x", type=float, default=-2.0)
ap.add_argument("--png")
a = ap.parse_args()
data = json.load(open(os.path.join(T, "data", "fig8_u_xm2.json")))

def section(r, x0):
    """the columns of the cross section x* = x0 (the mean of the two columns beside x0 if x0 is on a face): y*, the
    depth (m) and u3/U0 (ny, nz), sorted by y*"""
    main = r["y"] >= 0
    dist = np.abs(r["x"] - x0)
    sel = main & (dist < 0.5*r["dx"] + 1e-6)
    ys = np.unique(np.round(r["y"][sel], 6))
    Y, H, U = [], [], []
    for y in ys:
        k = sel & (np.abs(r["y"] - y) < 1e-5)
        Y.append(y); H.append(r["h"][k].mean()); U.append(r["u3"][k].mean(0)/U0)
    return np.array(Y), np.array(H), np.array(U)

def interp(Y, H, U, yq, zq):
    """u/U0 at (y/W, z/W): linear in z within the columns (layers of equal thickness, z from the bed), then in y"""
    nz = U.shape[1]
    def col(j):
        zc = (np.arange(nz) + 0.5)/nz*H[j]/W
        return np.interp(zq, zc, U[j])
    j = np.clip(np.searchsorted(Y, yq) - 1, 0, len(Y) - 2)
    w = np.clip((yq - Y[j])/(Y[j+1] - Y[j]), 0, 1)
    return (1 - w)*col(j) + w*col(j+1)

def width(ys, us, thr=-0.5):
    """first y with u <= thr (linear between the points)"""
    k = np.nonzero(us <= thr)[0]
    if len(k) == 0:
        return np.nan
    k = k[0]
    if k == 0:
        return ys[0]
    return ys[k-1] + (thr - us[k-1])*(ys[k] - ys[k-1])/(us[k] - us[k-1])

res = {}
for run, case in (("RUN5", a.run5), ("RUN2", a.run2)):
    if not case:
        continue
    r = load_any(case)
    Y, H, U = section(r, a.x)
    d = np.array(data[run])
    m = np.array([interp(Y, H, U, y, z) for y, z, _, _ in d])
    e = np.clip(m, -1.4, 0.1) - d[:, 2]
    zs = np.unique(d[:, 1])
    zb, zt = zs[1], zs[-1]                       # near the bed (z/W 0.03) and the top row of the measurements
    wd = [width(d[d[:, 1] == z][:, 0], d[d[:, 1] == z][:, 2]) for z in (zb, zt)]
    yy = np.linspace(Y[0], Y[-1], 400)
    wm = [width(yy, np.array([interp(Y, H, U, y, z) for y in yy])) for z in (zb, zt)]
    print("%s %-12s t %.0f s  nz %d  depth at x*=%g: %.4f m  rms (clipped) %.3f  mean %+.3f  |  u/U0 = -0.5 at y/W: "
          "bed (z/W %.2f) meas %.2f model %.2f, top (z/W %.2f) meas %.2f model %.2f  |  model min u/U0 %.2f, max %+.2f"
          % (run, os.path.basename(os.path.normpath(case)), r["simtime"], U.shape[1], a.x, H.mean(), np.sqrt((e**2).mean()),
             e.mean(), zb, wd[0], wm[0], zt, wd[1], wm[1], U.min(), U.max()))
    res[run] = (Y, H, U, d, zb, zt, wd, wm)

if a.png and res:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm
    C = ["#2B4C7E", "#4F8FCF", "#E07B39", "#3A9A5B", "#6B6B6B"]
    lev = np.arange(-1.45, 0.16, 0.1)
    cm = plt.get_cmap("RdYlBu", len(lev) - 1)
    norm = BoundaryNorm(lev, cm.N)
    fig, ax = plt.subplots(len(res), 3, figsize=(14, 3.4*len(res)), squeeze=False,
                           gridspec_kw={"width_ratios": [1, 1, 0.9]}, layout="constrained")
    for i, (run, (Y, H, U, d, zb, zt, wd, wm)) in enumerate(res.items()):
        ys, zs = np.unique(d[:, 0]), np.unique(d[:, 1])
        G = np.full((len(zs), len(ys)), np.nan)
        for y, z, u, _ in d:
            G[np.searchsorted(zs, z), np.searchsorted(ys, y)] = u
        ax[i, 0].pcolormesh(ys, zs, G, cmap=cm, norm=norm, shading="nearest")
        ax[i, 0].set_title("%s measured (Huang et al. 2002)" % run, fontsize=9)
        nz = U.shape[1]
        yg = np.linspace(0.005, 0.995, 100)
        zg = np.linspace(0.002, H.mean()/W - 0.002, 60)
        M = np.array([[interp(Y, H, U, y, z) for y in yg] for z in zg])
        pc = ax[i, 1].pcolormesh(yg, zg, np.clip(M, -1.4, 0.1), cmap=cm, norm=norm, shading="nearest")
        ax[i, 1].set_title("%s NHFLOW, %d layers" % (run, nz), fontsize=9)
        for k in (0, 1):
            ax[i, k].set_xlim(0, 1); ax[i, k].set_ylim(0, 0.36)
            ax[i, k].set_xlabel("y/W  (0: branch side)"); ax[i, k].set_ylabel("z/W")
            ax[i, k].set_aspect("equal")
        fig.colorbar(pc, ax=ax[i, 1], shrink=0.9, label="u/U0")
        yy = np.linspace(Y[0], Y[-1], 200)
        for c, z in enumerate((zb, 0.15, zt)):
            s = d[np.abs(d[:, 1] - z) < 1e-6]
            ax[i, 2].plot(s[:, 0], s[:, 2], "o", mfc="none", ms=3, color=C[c])
            ax[i, 2].plot(yy, [interp(Y, H, U, y, z) for y in yy], "-", lw=1.4, color=C[c], label="z/W = %.2f" % z)
        ax[i, 2].axhline(-1.4, color=C[4], lw=0.8, ls=":")
        ax[i, 2].set_xlim(0, 1); ax[i, 2].set_ylim(-1.8, 0.4)
        ax[i, 2].set_xlabel("y/W"); ax[i, 2].set_ylabel("u/U0")
        ax[i, 2].grid(alpha=0.3)
        for s in ("top", "right"):
            ax[i, 2].spines[s].set_visible(False)
        ax[i, 2].legend(fontsize=8, frameon=False, loc="upper right", title="lines: model, circles: measured",
                        title_fontsize=7)
    fig.suptitle("u/U0 at x/W = %g (positive upstream, U0 = %.3f m/s)" % (a.x, U0), fontsize=10)
    fig.savefig(a.png, dpi=150)
    print("wrote", a.png)
