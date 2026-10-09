# Animation of the FW-H observer planes of validation case 15 (circles) over the Curle dipole field of the
# box force (background): near field +-40 m (p' r^2) and acoustic field +-40 km (p' r), as MP4 (ffmpeg).
#   python sphere_fwh_animation.py OUTDIR [t0 [periods [St]]] [--label "sphere Re = 300, D/20"]
# in the run directory; default t0 = 150 s, 3 periods of St = 0.1357 (U = D = 1), frames every 0.1 s.
# Writes OUTDIR/fwh_planes.mp4 and fwh_planes_still.png. For PowerPoint/Keynote a looping GIF of one period:
#   ffmpeg -t <period*10/30> -i fwh_planes.mp4 -vf "fps=15,scale=1200:-1:flags=lanczos,split[a][b];
#     [a]palettegen=max_colors=128[p];[b][p]paletteuse=dither=none" -loop 0 fwh_planes.gif
import math, os, sys, collections
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
from matplotlib.patches import Rectangle, Circle

label = "sphere Re = 300"
argv = sys.argv[1:]
if "--label" in argv:
    i = argv.index("--label")
    label = argv[i+1]
    del argv[i:i+2]
OUT = argv[0] if argv else "slides"
os.makedirs(OUT, exist_ok=True)
t_start = float(argv[1]) if len(argv) > 1 else 150.0
periods = float(argv[2]) if len(argv) > 2 else 3.0
St = float(argv[3]) if len(argv) > 3 else 0.1357
for ff in ("/opt/homebrew/bin/ffmpeg", "/usr/local/bin/ffmpeg"):
    if os.path.exists(ff):
        plt.rcParams["animation.ffmpeg_path"] = ff
DARK, LIGHT, BODY, ORANGE = "#0E2131", "#F5F3EE", "#B8C7D3", "#E8894F"
plt.rcParams.update({"font.family": "Helvetica Neue", "font.size": 15, "text.color": LIGHT,
                     "axes.labelcolor": LIGHT, "xtick.color": BODY, "ytick.color": BODY,
                     "axes.edgecolor": BODY, "figure.facecolor": DARK, "savefig.facecolor": DARK})
c0 = 1500.0
obs = []
for l in open("ctrl.txt"):
    s = l.split()
    if len(s) >= 5 and s[0] == "U" and s[1] == "30":
        obs.append(np.array([float(x) for x in s[2:5]]))


def table(fn, ncol):
    d = {}
    for l in open(fn):
        s = l.split()
        if len(s) < ncol:
            continue
        try:
            v = [float(x) for x in s[:ncol]]
        except ValueError:
            continue
        d[v[0]] = v[1:]
    t = np.array(sorted(d))
    return t, np.array([d[x] for x in t])


def signal(n):
    d = {}
    for l in open(f"REEF3D_CFD_Acoustics/REEF3D-CFD-FWH-Observer-{n+1}.dat"):
        s = l.split()
        if len(s) != 2:
            continue
        try:
            d[float(s[0])] = float(s[1])
        except ValueError:
            pass
    t = np.array(sorted(d))
    return t, np.array([d[x] for x in t])


tb, B = table("REEF3D_CFD_Acoustics/REEF3D-CFD-FWH-Box.dat", 13)
FB = -(B[:, 0:3] + B[:, 3:6] + np.gradient(B[:, 9:12], tb, axis=0))
FB -= FB[tb > 0.5*tb[-1]].mean(axis=0)
dFB = np.gradient(FB, tb, axis=0)

T0, T1, DT = t_start, t_start + periods/St, 0.1
times = np.arange(T0, T1, DT)

planes = collections.defaultdict(dict)
for n, (x, y, z) in enumerate(obs):
    for L in (5.0, 5000.0):
        i, j = round(x/L), round(y/L)
        if abs(x-i*L) < 1e-6*L and abs(y-j*L) < 1e-6*L and abs(i) <= 8 and abs(j) <= 8 and abs(z) < 0.1 and max(abs(i), abs(j)) > 0:
            planes[L][(i, j)] = n


def curle_field(X, Y, t):
    R = np.hypot(X, Y)
    tau = t - R/c0
    F = [np.interp(tau, tb, FB[:, i]) for i in range(2)]
    dF = [np.interp(tau, tb, dFB[:, i]) for i in range(2)]
    return -((X*F[0] + Y*F[1])/R**3 + (X*dF[0] + Y*dF[1])/(c0*R**2))/(4*np.pi)


panels = []
for L, pw, unit, title in ((5.0, 2, "m", "hydrodynamic near field  ±40 m:  p′·r²"),
                           (5000.0, 1, "km", "acoustic far field  ±40 km:  p′·r")):
    E = 8*L
    s = np.linspace(-E, E, 320)
    X, Y = np.meshgrid(s, s)
    R = np.hypot(X, Y)
    keys = sorted(planes[L])
    xs = np.array([k[0]*L for k in keys])
    ys = np.array([k[1]*L for k in keys])
    rs = np.hypot(xs, ys)
    vals = []
    for k in keys:
        t, p = signal(planes[L][k])
        w = (t > T0 - 30.0) & (t < T1 + 30.0)
        p = p - p[w].mean()
        vals.append(np.interp(times, t, p))
    vals = np.array(vals)*rs[:, None]**pw
    vmax = max(np.abs(curle_field(X, Y, t)*R**pw)[R > 0.25*E].max() for t in times[::10])
    panels.append(dict(L=L, E=E, X=X, Y=Y, R=R, pw=pw, unit=unit, title=title, xs=xs, ys=ys, vals=vals, vmax=vmax))

fig, axs = plt.subplots(1, 2, figsize=(16, 8.0))
fig.subplots_adjust(left=0.05, right=0.97, top=0.84, bottom=0.08, wspace=0.18)
art = []
for a, P in zip(axs, panels):
    sc = 1e3 if P["unit"] == "km" else 1.0
    a.set_facecolor(DARK)
    im = a.pcolormesh(P["X"]/sc, P["Y"]/sc, curle_field(P["X"], P["Y"], T0)*P["R"]**P["pw"], cmap="RdBu_r",
                      vmin=-P["vmax"], vmax=P["vmax"], shading="auto", rasterized=True)
    pts = a.scatter(P["xs"]/sc, P["ys"]/sc, c=P["vals"][:, 0], cmap="RdBu_r", vmin=-P["vmax"], vmax=P["vmax"],
                    s=70 if P["unit"] == "m" else 80, edgecolors=DARK, linewidths=1.2, zorder=3)
    if P["unit"] == "m":
        a.add_patch(Rectangle((-1.5, -1.5), 4, 3, fill=False, ec=DARK, lw=1.5, zorder=4))
        a.add_patch(Circle((0, 0), 0.5, color=DARK, zorder=4))
        a.annotate("FW-H box + sphere", (0.5, -1.6), (12, -37), color=LIGHT, fontsize=14,
                   arrowprops=dict(arrowstyle="-", color=LIGHT, lw=1))
        a.text(-40, 41.5, "flow in +x", color=LIGHT, fontsize=15)
    else:
        a.add_patch(Circle((0, 0), c0/(2*np.pi*St)/1e3, fill=False, ec=LIGHT, lw=1, ls=":", zorder=4))
    a.set_aspect("equal")
    a.set_xlim(-P["E"]/sc, P["E"]/sc)
    a.set_ylim(-P["E"]/sc, P["E"]/sc)
    a.set_xlabel(f"x [{P['unit']}]")
    a.set_ylabel(f"y [{P['unit']}]")
    a.set_title(P["title"], color=LIGHT, fontsize=17, pad=12)
    for sp in a.spines.values():
        sp.set_visible(False)
    art.append((im, pts, P))
head = fig.text(0.05, 0.93, "", fontsize=19, color=LIGHT)
fig.text(0.97, 0.93, "background: Curle dipole of the box force   ·   circles: FW-H observers", fontsize=15, color=BODY, ha="right")

writer = FFMpegWriter(fps=30, codec="libx264", bitrate=6000, extra_args=["-pix_fmt", "yuv420p", "-preset", "slow"])
with writer.saving(fig, os.path.join(OUT, "fwh_planes.mp4"), dpi=100):
    for f, t in enumerate(times):
        for im, pts, P in art:
            im.set_array((curle_field(P["X"], P["Y"], t)*P["R"]**P["pw"]).ravel())
            pts.set_array(P["vals"][:, f])
        head.set_text(f"{label} — t = {t:6.1f} s")
        writer.grab_frame()
        if f == 0:
            fig.savefig(os.path.join(OUT, "fwh_planes_still.png"), dpi=100)
print(len(times), "frames")
