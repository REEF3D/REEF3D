# The Weber junction runs of weber_case.py against the measurements: the water surface along the main channel
# (Huang et al. 2002, Fig. 5, RUN5), its minimum in the separation zone and its rise across the junction (0.95 and
# +5.0 % for RUN5, 0.92 and +8.8 % for RUN2), the separation zone (reversed depth-averaged u) and the surface map.
# The state is the AMR dump of the run (REEF3D_AMR_DUMP=<case>/amr); z/H0 is the depth over the mean depth across
# x* = -4 (the first points of Fig. 5 are 1.000 there).
# usage: python3 weber_compare.py case [case ...] [--np 1] [--png out.png] [--labels a,b]
import argparse, glob, json, os, re, sys, numpy as np
T = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(T, "../sediment_curvilinear"))
from p5lib import amr_grids2

ap = argparse.ArgumentParser()
ap.add_argument("cases", nargs="+")
ap.add_argument("--png")
ap.add_argument("--labels")
a = ap.parse_args()
FT = 0.3048
W, Lm, Lb, Lu = 3*FT, 72*FT, 12*FT, 18*FT

def read_blocks(path):
    """blocks.dat: the node arrays (ni + 1, nj + 1) of x and y per block"""
    t = open(path).read().split()
    k, out = 1, []
    for _ in range(int(t[0])):
        ni, nj = int(t[k]), int(t[k+1]); k += 2
        a = np.array(t[k:k+3*(ni+1)*(nj+1)], float).reshape(nj+1, ni+1, 3); k += 3*(ni+1)*(nj+1)
        out.append((a[:, :, 0].T, a[:, :, 1].T))
    return out

def load(d):
    """cell centres (x*, y* of Huang et al.), eta, depth and the velocity along x*, y* of all level-0 cells (the AMR
    dump has the Cartesian velocities); the cells from blocks.dat, block 0's from the native grid less its ring"""
    nodes = read_blocks(os.path.join(d, "blocks.dat"))
    ring = [0, 0, 0, 0]
    for l in open(os.path.join(d, "reef.log"), errors="replace"):
        m = re.search(r"ring (\d+) (\d+) (\d+) (\d+)", l)
        if m:
            ring = [int(v) for v in m.groups()]
    npr = len(glob.glob(os.path.join(d, "amr", "amr_r*.bin")))
    X, Y, Hd, U, V = [], [], [], [], []
    for g in amr_grids2(d, npr):
        if g["lev"] != 0:
            continue
        # block 0: its cells from (I0, J0) of the native grid; the blocks 1..: their own cells (I0, J0 place the level-0
        # grid of a block in the index space of the hierarchy)
        o = (g["I0"] - ring[0], g["J0"] - ring[2]) if g["blk"] == 0 else (0, 0)
        nx_, ny_ = nodes[g["blk"]][0].shape[0]-1, nodes[g["blk"]][0].shape[1]-1
        I = o[0] + np.arange(g["nx"])[:, None] + 0*np.arange(g["ny"])[None, :]
        J = o[1] + 0*np.arange(g["nx"])[:, None] + np.arange(g["ny"])[None, :]
        keep = (I >= 0) & (I < nx_) & (J >= 0) & (J < ny_)
        x, y = nodes[g["blk"]]
        cx = 0.25*(x[:-1, :-1] + x[1:, :-1] + x[:-1, 1:] + x[1:, 1:])
        cy = 0.25*(y[:-1, :-1] + y[1:, :-1] + y[:-1, 1:] + y[1:, 1:])
        X.append(cx[I[keep], J[keep]]); Y.append(cy[I[keep], J[keep]])
        Hd.append(g["WL"][keep]); U.append(g["U"][keep]); V.append(g["V"][keep])
    X, Y, Hd, U, V = (np.concatenate(q) for q in (X, Y, Hd, U, V))
    dx = np.hypot(*(np.diff(nodes[0][0][:2, 0]), np.diff(nodes[0][1][:2, 0])))[0]
    xs, ys = (Lu - X)/W, (W - Y)/W                 # Huang's x*, y*
    return dict(x=xs, y=ys, h=Hd, u=-U, v=-V, dx=dx/W, simtime=g["simtime"])

def read_vtp(fn):
    """the points and the point data of a VTK PolyData file of REEF3D (raw appended data, UInt32 headers)"""
    buf = open(fn, "rb").read()
    k = buf.index(b"<AppendedData")
    head = buf[:k].decode(errors="replace")
    start = buf.index(b"_", k) + 1
    out = {}
    for m in re.finditer(r'<DataArray type="(\w+)"(?: Name="(\w+)")?(?: NumberOfComponents="(\d)")? format="appended" offset="(\d+)"', head):
        typ, name, nc, off = m.group(1), m.group(2) or "Points", int(m.group(3) or 1), int(m.group(4))
        n = int(np.frombuffer(buf, "<u4", 1, start+off)[0])
        dt = {"Float32": "<f4", "Float64": "<f8", "Int32": "<i4"}[typ]
        out[name] = np.frombuffer(buf, dt, n//np.dtype(dt).itemsize, start+off+4).reshape(-1, nc) if nc > 1 else \
            np.frombuffer(buf, dt, n//np.dtype(dt).itemsize, start+off+4)
    return out

def load_cart(d):
    """the Cartesian reference of weber_cart.py: the last VTP of the free surface (point data), the points inside
    the channels"""
    meta = json.load(open(os.path.join(d, "case.json")))
    last = sorted(glob.glob(os.path.join(d, "REEF3D_SFLOW_VTP_FSF", "*.pvtp")))[-1]
    pieces = re.findall(r'Source="([^"]+)"', open(last).read())
    X, Y, Hd, U, V, t = [], [], [], [], [], 0.0
    for pc in pieces:
        v = read_vtp(os.path.join(os.path.dirname(last), pc))
        X.append(v["Points"][:, 0]); Y.append(v["Points"][:, 1]); Hd.append(v["waterlevel"])
        U.append(v["velocity"][:, 0]); V.append(v["velocity"][:, 1])
    X, Y, Hd, U, V = (np.concatenate(q).astype(float) for q in (X, Y, Hd, U, V))
    t = float(re.search(rb'Name="TimeValue"[^>]*>\s*([0-9.eE+-]+)', open(os.path.join(os.path.dirname(last), pieces[0]), "rb").read()).group(1))
    e = 1e-6
    inside = (Y <= W+e) | ((X >= Lu-e) & (X <= Lu+W+e))
    X, Y, Hd, U, V = (q[inside] for q in (X, Y, Hd, U, V))
    xs, ys = (Lu - X)/W, (W - Y)/W
    return dict(x=xs, y=ys, h=Hd, u=-U, v=-V, dx=meta["dx"]/W, simtime=t, note=meta["note"])

def at(r, xs, ys):
    """nearest cell"""
    k = np.argmin((r["x"]-xs)**2 + (r["y"]-ys)**2)
    return k

exp = json.load(open(os.path.join(T, "data", "run5_surface_fig5.json")))["profiles"]
R = [load_cart(d) if os.path.exists(os.path.join(d, "case.json")) else load(d) for d in a.cases]
labels = a.labels.split(",") if a.labels else [os.path.basename(os.path.normpath(d)) for d in a.cases]
for d, r, lab in zip(a.cases, R, labels):
    main = r["y"] >= 0
    sec = lambda x0: main & (np.abs(r["x"]-x0) < 0.5*r["dx"]+1e-9)
    H0 = r["h"][sec(-4.0 + 0.5*r["dx"]) | sec(-4.0 - 0.5*r["dx"])].mean()
    z = r["h"]/H0
    up = z[sec(1.0 + 0.5*r["dx"]) | sec(1.0 - 0.5*r["dx"])].mean()
    junc = main & (r["x"] < -1.15) & (r["x"] > -4)     # below the corner cell of the branch
    kmin = np.argmin(np.where(junc, z, 9))
    # separation zone: reversed depth-averaged u (u* > 0: towards x* > 0, upstream) in the main channel below x* = -1
    rev = main & (r["x"] < -1) & (r["u"] > 0)
    Ls = (-1 - r["x"][rev].min()) if rev.any() else 0.0
    Ws = r["y"][rev].max() + 0.5*r["dx"] if rev.any() else 0.0
    h7 = r["h"][sec(-7.0 + 0.5*r["dx"]) | sec(-7.0 - 0.5*r["dx"])].mean()
    print("%-10s t %.0f s  H0(x*=-4) %.4f m  h(x*=-7) %.4f  upstream (x*=1) %.4f (+%.1f %%)  min %.4f at x* %.2f y* %.2f  "
          "reversed flow: length %.2f W, width %.2f W" % (lab, r["simtime"], H0, h7, up, 100*(up-1), z[kmin], r["x"][kmin],
                                                       r["y"][kmin], Ls, Ws))
    r["z"], r["H0"] = z, H0
    if "q* = 0.75" in (r.get("note") or json.load(open(os.path.join(d, "blocks.json")))["note"]):
        err = []
        for k, pts in exp.items():
            yv = float(k)
            for xv, zv in pts:
                err.append(z[at(r, xv, yv)] - zv)
        err = np.array(err)
        print("%-10s Fig. 5 surface: %d points, mean %+.4f, rms %.4f, max |%.4f|" % ("", len(err), err.mean(),
                                                                                  np.sqrt((err**2).mean()), np.abs(err).max()))

if a.png:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    C = ["#2B4C7E", "#4F8FCF", "#E07B39", "#3A9A5B", "#6B6B6B"]
    fig, ax = plt.subplots(len(exp), 1, figsize=(7.5, 10), sharex=True)
    for i, (k, pts) in enumerate(exp.items()):
        yv = float(k)
        p = np.array(pts)
        ax[i].plot(p[:, 0], p[:, 1], "o", mfc="none", color=C[4], ms=4, label="measured (Fig. 5)")
        for m, (r, lab) in enumerate(zip(R, labels)):
            line = (np.abs(r["y"] - yv) < 0.5*r["dx"]+1e-9) & (r["x"] > -4.5) & (r["x"] < 1.5)
            if not line.any():
                yy = r["y"][np.argmin(np.abs(r["y"]-yv))]
                line = (r["y"] == yy) & (r["x"] > -4.5) & (r["x"] < 1.5)
            o = np.argsort(r["x"][line])
            ax[i].plot(r["x"][line][o], r["z"][line][o], "-", color=C[m % 4], lw=1.4, label=lab)
        ax[i].set_ylim(0.92, 1.10)
        ax[i].text(0.98, 0.1, "y/W = %s" % k, transform=ax[i].transAxes, ha="right")
        ax[i].grid(alpha=0.3)
        for s in ("top", "right"):
            ax[i].spines[s].set_visible(False)
        ax[i].axvspan(-1, 0, color="#dddddd", alpha=0.4, lw=0)
    ax[0].legend(loc="upper left", fontsize=8, frameon=False)
    ax[-1].set_xlabel("x/W  (positive upstream, branch at -1 .. 0)")
    ax[len(exp)//2].set_ylabel("z / H0")
    ax[-1].set_xlim(-4.2, 1.2)
    fig.tight_layout()
    fig.savefig(a.png, dpi=150)
    print("wrote", a.png)
