# Animations of a KP505 open-water run from its LAGOON store (P 18), for slides:
#   python kp505_lagoon_anim.py [t0 t1] [--out DIR] [--fps 12]      in the run directory
# 1 kp505_meridional.mp4: plane z = 0, axial velocity and azimuthal vorticity, the blade sections
# 2 kp505_rotorplane.mp4: plane just upstream of the rotor (x = -0.04), pressure fluctuation p' =
#   p - mean, the rotating BPF pattern of the near field, with the probe positions (P 64)
# each with a still (PNG) and a looping GIF of one revolution (for PowerPoint/Keynote).
# Needs numpy, matplotlib, ffmpeg; lagoon_read.py next to this script.
import argparse, os, shutil, subprocess, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lagoon_read import lagoon

ap = argparse.ArgumentParser()
ap.add_argument("t", nargs="*", type=float)
ap.add_argument("--out", default="anim")
ap.add_argument("--fps", type=float, default=12.0)
args = ap.parse_args()
t0, t1 = (args.t + [0.4, 0.8])[:2] if len(args.t) < 2 else args.t[:2]
os.makedirs(args.out, exist_ok=True)

DARK, LIGHT, BODY = "#0E2131", "#F5F3EE", "#B8C7D3"
plt.rcParams.update({"font.family": "Helvetica Neue", "font.size": 14, "text.color": LIGHT, "axes.labelcolor": LIGHT,
                     "xtick.color": BODY, "ytick.color": BODY, "axes.edgecolor": BODY,
                     "figure.facecolor": DARK, "savefig.facecolor": DARK})

n_rps = U = None
probes = []
for l in open("ctrl.txt"):
    s = l.split()
    if len(s) >= 3 and s[0] == "X" and s[1] == "211": n_rps = abs(float(s[2]))/(2*np.pi)
    if len(s) >= 3 and s[0] == "U" and s[1] == "50": U = float(s[2])
    if len(s) >= 5 and s[0] == "P" and s[1] == "64": probes.append([float(v) for v in s[2:5]])
probes = np.array(probes)

S = lagoon("REEF3D_CFD.lagoon")
x, y, z = S.coords()
frames = [k for k, t in enumerate(S.times) if t0 - 1e-9 <= t <= t1 + 1e-9]
print(f"{len(frames)} frames {S.times[frames[0]]:.3f}-{S.times[frames[-1]]:.3f} s, n = {n_rps} 1/s, U = {U} m/s")

# meridional window and rotor-plane window
ix = (x >= -0.15) & (x <= 0.5)
iy = (y >= -0.25) & (y <= 0.25)
k0 = int(np.argmin(abs(z)))
ixp = int(np.argmin(abs(x + 0.04)))
iyz = (y >= -0.4) & (y <= 0.4)


def section(V, T, zc=0.0):
    # segments of the triangles cut by the plane z = zc
    P = V[T]                                   # (n, 3, 3)
    d = P[:, :, 2] - zc
    seg = []
    for a, b in ((0, 1), (1, 2), (2, 0)):
        m = d[:, a]*d[:, b] < 0
        w = d[m, a]/(d[m, a] - d[m, b])
        seg.append((m, P[m, a] + w[:, None]*(P[m, b] - P[m, a])))
    out = []
    pts = {}
    for m, q in seg:
        idx = np.where(m)[0]
        for i, pt in zip(idx, q):
            pts.setdefault(i, []).append(pt[:2])
    for v in pts.values():
        if len(v) == 2:
            out.append(v)
    return np.array(out)


# mean pressure in the rotor plane (for p')
pm = 0.0
for k in frames:
    pm = pm + S.volume("pressure", k)[:, :, ixp]
pm /= len(frames)


def render(kind):
    fn = os.path.join(args.out, f"kp505_{kind}")
    tmp = fn + "_frames"
    os.makedirs(tmp, exist_ok=True)
    for j, k in enumerate(frames):
        t = S.times[k]
        V, T = S.body(k)
        if kind == "meridional":
            u = S.volume("velocity", k)
            ux = u[k0, :, :, 0][np.ix_(iy, ix)]
            # azimuthal vorticity in z = 0: dv/dx - du/dy
            uu = u[k0, :, :, 0]; vv = u[k0, :, :, 1]
            om = np.gradient(vv, x, axis=1) - np.gradient(uu, y, axis=0)
            om = om[np.ix_(iy, ix)]
            fig, ax = plt.subplots(2, 1, figsize=(10.5, 8.6), sharex=True)
            fig.subplots_adjust(left=0.1, right=0.95, top=0.91, bottom=0.08, hspace=0.12)
            im0 = ax[0].pcolormesh(x[ix], y[iy], ux, cmap="RdBu_r", vmin=0.0, vmax=2*U, shading="auto")
            im1 = ax[1].pcolormesh(x[ix], y[iy], om, cmap="PuOr_r", vmin=-150, vmax=150, shading="auto")
            seg = section(V, T)
            for a, im, lab in ((ax[0], im0, "axial velocity u [m/s]"), (ax[1], im1, "azimuthal vorticity [1/s]")):
                if len(seg):
                    a.add_collection(LineCollection(seg, colors="k", linewidths=1.6))
                a.set_aspect("equal")
                a.set_ylabel("y [m]")
                cb = fig.colorbar(im, ax=a, shrink=0.85, pad=0.01)
                cb.set_label(lab)
            ax[1].set_xlabel("x [m]   (flow in +x)")
            fig.text(0.1, 0.95, f"KP505, J = {U/(n_rps*0.25):.2f}, plane z = 0      t = {t:.3f} s", fontsize=17)
        else:
            p = S.volume("pressure", k)[:, :, ixp] - pm
            sub = p[np.ix_(iyz, iyz)]
            fig, ax = plt.subplots(figsize=(9.2, 8))
            fig.subplots_adjust(left=0.1, right=0.9, top=0.88, bottom=0.08)
            vm = 600.0
            im = ax.pcolormesh(y[iyz], z[iyz], sub, cmap="RdBu_r", vmin=-vm, vmax=vm, shading="auto")
            # blade outline projected on the plane (vertices upstream of the hub, every 20th)
            ax.plot(V[::5, 1], V[::5, 2], ",", color="k", alpha=0.3)
            ax.add_patch(plt.Circle((0, 0), 0.125, fill=False, ec=LIGHT, ls=":", lw=1))
            for q, pr in enumerate(probes):
                if abs(pr[0] - 0.003) < 0.01:
                    ax.plot(pr[1], pr[2], "o", ms=9, mfc="none", mec="#E8894F", mew=2)
            ax.set_aspect("equal")
            ax.set_xlabel("y [m]")
            ax.set_ylabel("z [m]")
            cb = fig.colorbar(im, ax=ax, shrink=0.85)
            cb.set_label("p' [Pa]")
            fig.text(0.1, 0.92, f"x = {x[ixp]:.3f} m (8 mm upstream of the blades): p' = p - mean      t = {t:.3f} s", fontsize=15)
        fig.savefig(os.path.join(tmp, f"f{j:04d}.png"), dpi=100)
        if j == 0:
            fig.savefig(fn + "_still.png", dpi=130)
        plt.close(fig)
        print(kind, j+1, "/", len(frames), flush=True)
    ff = shutil.which("ffmpeg") or "/opt/homebrew/bin/ffmpeg"
    subprocess.run([ff, "-loglevel", "error", "-y", "-framerate", f"{args.fps:g}", "-i", os.path.join(tmp, "f%04d.png"),
                    "-c:v", "libx264", "-pix_fmt", "yuv420p", "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", fn + ".mp4"], check=True)
    # GIF of one revolution
    nrev = int(round(1.0/(n_rps*(S.times[frames[1]] - S.times[frames[0]]))))
    subprocess.run([ff, "-loglevel", "error", "-y", "-framerate", f"{args.fps:g}", "-i", os.path.join(tmp, "f%04d.png"),
                    "-frames:v", str(max(nrev, 2)), "-vf",
                    "scale=900:-1:flags=lanczos,split[a][b];[a]palettegen=max_colors=128[p];[b][p]paletteuse=dither=none",
                    "-loop", "0", fn + ".gif"], check=True)
    print("written", fn + ".mp4", fn + ".gif", fn + "_still.png")


render("meridional")
render("rotorplane")
