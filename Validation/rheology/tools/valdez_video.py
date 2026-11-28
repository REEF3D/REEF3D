#!/usr/bin/env python3
# valdez_video.py: video of a Port Valdez run (case 36, NHFLOW or SFLOW free-surface output with the slide layer):
# water level over the bay and at the head of the bay, slide thickness, gauges at the waterfront.
#
#   python3 valdez_video.py <run> <case dir with reference.npz> <out.mp4> [--fps 12] [--label text]
#
# Needs matplotlib and ffmpeg. Frames are written to <out>_frames/ and joined with ffmpeg (H.264, yuv420p).
import argparse, glob, os, re, subprocess, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "turbulence", "tools"))
sys.path.insert(0, os.path.dirname(__file__))
from r3read import read_vtk


def frames(run):
    for d, pat in (("REEF3D_NHFLOW_VTP_FSF", r"FSF-(\d+)-(\d+)\.vtp$"), ("REEF3D_SFLOW_VTP_FSF", r"FSF-(\d+)-(\d+)\.vtp$")):
        fs = glob.glob(os.path.join(run, d, "*.vtp"))
        if fs:
            byn = {}
            for f in fs:
                m = re.search(pat, f)
                if m:
                    byn.setdefault(int(m.group(1)), []).append(f)
            return d.startswith("REEF3D_NHFLOW"), [byn[k] for k in sorted(byn)]
    raise SystemExit("no free-surface output in %s" % run)


def grid(pieces, nx, ny, dx, keys):
    """node values on the (nx+1, ny+1) node grid"""
    out = {k: np.full((ny + 1, nx + 1), np.nan) for k in keys}
    t = None
    for f in pieces:
        d = read_vtk(f)
        t = d.get("time", t)
        P = d["points"]
        i = np.clip(np.rint(P[:, 0]/dx).astype(int), 0, nx); j = np.clip(np.rint(P[:, 1]/dx).astype(int), 0, ny)
        for k in keys:
            if k in d:
                v = d[k]
                out[k][j, i] = v if v.ndim == 1 else np.linalg.norm(v, axis=1)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run"); ap.add_argument("case"); ap.add_argument("out")
    ap.add_argument("--fps", type=int, default=12)
    ap.add_argument("--label", default="")
    ap.add_argument("--keep", action="store_true", help="keep frames already rendered in <out>_frames")
    ap.add_argument("--dt", type=float, default=None, help="time between frames (s); default from P 182 in ctrl.txt")
    a = ap.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LightSource, LinearSegmentedColormap
    # water level: deep blue - pale blue (still water) - red
    wcm = LinearSegmentedColormap.from_list("wl", ["#08306B", "#4F8FCF", "#D6E9F8", "#F4A582", "#B2182B"])
    from rheo_eval import _gauges
    R = np.load(os.path.join(a.case, "reference.npz"))
    dx, bed = float(R["dx"]), R["bed"]
    ny, nx = bed.shape
    xn, yn = dx*np.arange(nx + 1), dx*np.arange(ny + 1)
    # terrain at the nodes (the mean of the cells around)
    bp = np.pad(bed, 1, mode="edge")
    bedn = 0.25*(bp[:-1, :-1] + bp[1:, :-1] + bp[:-1, 1:] + bp[1:, 1:])
    hs = LightSource(azdeg=315, altdeg=45).hillshade(np.maximum(bedn, 0.0), vert_exag=2, dx=dx, dy=dx)
    nh, fr = frames(a.run)
    dt = a.dt
    if dt is None:
        dt = 30.0
        for ln in open(os.path.join(a.run, "ctrl.txt")):
            w = ln.split()
            if len(w) >= 3 and w[0] == "P" and w[1] == "182":
                dt = float(w[2])
    t_g, G = _gauges(a.run)
    tend = None
    for ln in open(os.path.join(a.run, "ctrl.txt")):
        w = ln.split()
        if len(w) >= 3 and w[0] == "N" and w[1] == "41":
            tend = float(w[2])
    lines = [R[f] for f in R.files if f.startswith("line_InundationLine")]
    mck = [R[f] for f in R.files if f.startswith("line_McKinley")]
    fdir = a.out[:-4] + "_frames"
    os.makedirs(fdir, exist_ok=True)
    ext = [-0.5*dx, (nx + 0.5)*dx, -0.5*dx, (ny + 0.5)*dx]
    zoom = ((22300, 25800), (3200, 6900))
    for n, pcs in enumerate(fr):
        t = n*dt
        if a.keep and os.path.exists(os.path.join(fdir, "f%04d.png" % n)) and n < len(fr) - 1:
            continue
        g = grid(pcs, nx, ny, dx, ("eta", "WL", "slide_h", "velocity"))
        wet = np.isfinite(g["WL"]) & (g["WL"] > 0.05)
        eta = np.where(wet, g["eta"], np.nan)
        sh = np.where(np.isfinite(g["slide_h"]) & (g["slide_h"] > 0.5), g["slide_h"], np.nan)
        land = wet & (bedn > 0)
        fig = plt.figure(figsize=(16, 9), dpi=100)
        gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], width_ratios=[1.0, 1.35], hspace=0.25, wspace=0.12,
                              left=0.04, right=0.97, top=0.93, bottom=0.07)
        ax = fig.add_subplot(gs[0, :])
        ax.imshow(hs, origin="lower", extent=ext, cmap="gray", vmin=-0.3, vmax=1.3)
        im = ax.imshow(eta, origin="lower", extent=ext, cmap=wcm, vmin=-10, vmax=10)
        ims = ax.imshow(sh, origin="lower", extent=ext, cmap="copper_r", vmin=0, vmax=120, alpha=0.75)
        ax.contour(xn, yn, bedn, [0], colors="k", linewidths=0.5)
        for o in R["obs"]:
            ax.plot(o[0], o[1], "^", color="yellow", ms=6, mec="k")
        ax.add_patch(plt.Rectangle((zoom[0][0], zoom[1][0]), zoom[0][1] - zoom[0][0], zoom[1][1] - zoom[1][0], fill=False,
                                   ec="k", lw=0.8, ls="--"))
        ax.set_xlim(0, nx*dx); ax.set_ylim(0, ny*dx); ax.set_aspect("equal"); ax.tick_params(labelsize=8)
        ax.set_title("Port Valdez, 27 March 1964 (NTHMP benchmark 7)%s   t = %4.0f s (%.1f min)"
                     % ((" - " + a.label) if a.label else "", t, t/60), fontsize=12)
        cb = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.005, fraction=0.03); cb.set_label("water level (m)", fontsize=8)
        cb2 = fig.colorbar(ims, ax=ax, shrink=0.8, pad=0.06, fraction=0.03); cb2.set_label("slide thickness (m)", fontsize=8)
        cb.ax.tick_params(labelsize=7); cb2.ax.tick_params(labelsize=7)
        ax = fig.add_subplot(gs[1, 0])
        ax.imshow(hs, origin="lower", extent=ext, cmap="gray", vmin=-0.3, vmax=1.3)
        ax.imshow(np.where(~land, eta, np.nan), origin="lower", extent=ext, cmap=wcm, vmin=-10, vmax=10)
        dep = np.where(land, g["WL"], np.nan)
        im3 = ax.imshow(dep, origin="lower", extent=ext, cmap="magma_r", vmin=0, vmax=4)
        ax.imshow(sh, origin="lower", extent=ext, cmap="copper_r", vmin=0, vmax=120, alpha=0.6)
        ax.contour(xn, yn, bedn, [0], colors="k", linewidths=0.5)
        for q in lines:
            ax.plot(q[:, 0], q[:, 1], "-", color="gold", lw=1.5)
        for q in mck:
            ax.plot(q[:, 0], q[:, 1], "-", color="k", lw=1.2)
        ax.plot(R["obs"][0][0], R["obs"][0][1], "^", color="yellow", ms=7, mec="k")
        ax.set_xlim(*zoom[0]); ax.set_ylim(*zoom[1]); ax.set_aspect("equal"); ax.tick_params(labelsize=7)
        ax.set_title("head of the bay (old Valdez): flow depth on land; gold: observed inundation, black: McKinley St",
                     fontsize=8)
        cb3 = fig.colorbar(im3, ax=ax, shrink=0.8, pad=0.01); cb3.set_label("flow depth on land (m)", fontsize=8)
        ax = fig.add_subplot(gs[1, 1])
        if t_g is not None:
            nm = ["Valdez Hotel: flow depth x 5", "harbour", "Narrows light", "off Anderson Bay", "off Jackson Point", "mid-bay"]
            col = ["#E07B39", "#2B4C7E", "#3A9A5B", "#6B6B6B", "#8E44AD", "#4F8FCF"]
            o = R["obs"][0]
            zh = float(bed[int(o[1]//dx), int(o[0]//dx)])
            for k in range(min(G.shape[1], 6)):
                v = 5.0*np.maximum(G[:, k] - zh, 0.0) if k == 0 else G[:, k]
                ax.plot(t_g/60, v, "-", color=col[k], lw=1.0, alpha=0.25)
                m = t_g <= t
                ax.plot(t_g[m]/60, v[m], "-", color=col[k], lw=1.4 if k < 2 else 1.0, label=nm[k])
            ax.axvline(t/60, color="k", lw=0.8)
            ax.set_xlim(0, (tend or t_g.max())/60); ax.set_ylim(-15, 15)
            ax.axvspan(0, 5, color="k", alpha=0.05); ax.axvspan(10, 15, color="r", alpha=0.06)
            ax.text(0.15, 13.2, "1st wave at the waterfront\n(< 5 min, eyewitnesses)", fontsize=7, va="top")
            ax.text(10.15, 13.2, "2nd wave\n(10-15 min later)", fontsize=7, va="top")
            ax.legend(fontsize=7, loc="lower right", ncol=2)
            ax.set_xlabel("time (min)", fontsize=8); ax.set_ylabel("water level (m)", fontsize=8); ax.tick_params(labelsize=7)
            ax.grid(color="#DDDDDD", lw=0.6)
        fig.savefig(os.path.join(fdir, "f%04d.png" % n))
        plt.close(fig)
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(a.fps), "-i", os.path.join(fdir, "f%04d.png"),
                    "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "23",
                    a.out], check=True)
    print("%d frames, %s" % (len(fr), a.out))


if __name__ == "__main__":
    main()
