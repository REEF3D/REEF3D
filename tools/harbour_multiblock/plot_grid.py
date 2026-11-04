# P5d: the blocks of nidelva_harbour_grid.py over the water depth (schematic) and the river corridor in detail
# usage: python3 plot_grid.py griddir out.png
import sys, os, json, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLUE, LBLUE, ORANGE, GREEN, GREY = "#2B4C7E", "#4F8FCF", "#E07B39", "#3A9A5B", "#6B6B6B"
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
d, out = sys.argv[1], sys.argv[2]
g = np.load(os.path.join(d, "grid.npz"))
meta = json.load(open(os.path.join(d, "grid.json")))
WL = meta["WL"]
xs, ys, W, bed = g["xs"], g["ys"], g["W"], g["bed"]
dep = np.where(W, WL-bed, np.nan)

fig, axs = plt.subplots(1, 2, figsize=(12, 7), gridspec_kw={"width_ratios": [1.15, 1]})
ax = axs[0]
im = ax.imshow(dep.T, origin="lower", extent=[xs[0], xs[-1], ys[0], ys[-1]], cmap="Blues", vmin=0, vmax=25)
fig.colorbar(im, ax=ax, shrink=0.6, label="water depth (m, schematic)")
CX, CY, BX, BY = g["CX"], g["CY"], g["BX"], g["BY"]
def outline(X, Y, c, lw=1.2):
    ax.plot(np.r_[X[:, 0], X[-1, :], X[::-1, -1], X[0, ::-1]], np.r_[Y[:, 0], Y[-1, :], Y[::-1, -1], Y[0, ::-1]], c=c, lw=lw)
outline(CX, CY, ORANGE)
nb = meta.get("harbour_blocks", 1)
cuts = np.linspace(0, BX.shape[1]-1, nb+1).round().astype(int)
for b, (c0, c1) in enumerate(zip(cuts[:-1], cuts[1:])):
    outline(BX[:, c0:c1+1], BY[:, c0:c1+1], GREEN)
    ax.text(BX[-1, (c0+c1)//2], BY[-1, 0]+40, "block %d: water level (G 54 %d 2)" % (b+1, b+1), ha="center", color=GREEN, fontsize=8)
ax.annotate("block 0 (native): inflow Q (B 60 1, W 10)", xy=(CX[0, 5], CY[0, 5]), xytext=(1150, -2150), color=ORANGE,
            fontsize=8, arrowprops=dict(arrowstyle="->", color=ORANGE))
ax.annotate("interface", xy=(1230, -300), xytext=(1500, -650), fontsize=8, color=GREY, arrowprops=dict(arrowstyle="->", color=GREY))
ax.set_xlim(-600, 2500); ax.set_ylim(-2300, 2200); ax.set_aspect("equal")
ax.set_xlabel("x (m, local frame of georef.json)"); ax.set_ylabel("y (m)")
ax.set_title("Nidelva - Trondheim harbour: %d blocks, %d + %d cells" % (1+nb, CX.size - CX.shape[0] - CX.shape[1] + 1,
                                                                        (BX.shape[0]-1)*(BX.shape[1]-1)))
ax = axs[1]
ax.imshow(np.where(W, 1.0, np.nan).T, origin="lower", extent=[xs[0], xs[-1], ys[0], ys[-1]], cmap="Blues", vmin=0, vmax=3)
for k in range(CX.shape[0]):
    ax.plot(CX[k], CY[k], c=ORANGE, lw=0.4)
for k in range(CX.shape[1]):
    ax.plot(CX[:, k], CY[:, k], c=ORANGE, lw=0.4)
for k in range(0, BX.shape[0], 1):
    ax.plot(BX[k], BY[k], c=GREEN, lw=0.3)
for k in range(0, BX.shape[1], 1):
    ax.plot(BX[:, k], BY[:, k], c=GREEN, lw=0.3)
ax.set_xlim(900, 1500); ax.set_ylim(-950, -100); ax.set_aspect("equal")
ax.set_title("the corridor ends on the nodes of block 1 (10 m)")
fig.tight_layout()
fig.savefig(out, dpi=150)
