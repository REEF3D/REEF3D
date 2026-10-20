# Validation case cavitation 02 (hemispherical head form, Rouse & McNown 1948): surface pressure
# coefficient Cp(s/D) and cavity length from the VTU output, averaged over time and azimuth.
#
#   python headform_analysis.py [run dir ...] [--t0 T] [--t1 T] [--exp rouse_mcnown.csv] [--plot NAME]
#
# For each run directory (default: the current one): the points in two bands outside the body
# surface (distance nose: sqrt(x^2+r^2) - D/2, cylinder: r - D/2), 0.5-1.5 dx and 1.5-2.5 dx,
# binned in the surface coordinate from the nose tip, s = D/2 theta on the nose and s = pi/4 D + x
# on the cylinder; the wall pressure is extrapolated linearly from the two bands (the normal
# pressure gradient rho u^2/R at the curved nose is large: 0.05 D off the wall Cp is ~0.3 higher). Cp = (p - p_ref)/q
# with p_ref the mean of pressure probe 1 (upstream, P 64) over the window, q = rho U^2/2, U from
# W 10 over the inlet area; sigma = (p_ref - p_v)/q. Cavity: alpha_v > 0.1 in the band, contiguous
# from the first cavitating bin. --exp: csv with columns sigma, s/D, Cp (measured), drawn for the
# sigma closest to each run. Default window: the VTU frames of the last 5 s.
# Runs with mesh refinement (make_case.py --amr, P 18 1): the cell values of the finest patches from
# the LAGOON store REEF3D_CFD.lagoon (the bands in the refined box, at the patch spacing dx/2);
# needs lagoon_read.py (Regression/validation/fwh/scripts, or next to this script).
import argparse, glob, math, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(1, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "fwh", "scripts"))
from reef3d_vtu import steps, read_step


def store_frames(run):
    # AMR runs: (output number, time) of the AMR set and a reader of the finest patches of an output
    from lagoon_read import lagoon
    s = lagoon(f"{run}/REEF3D_CFD.lagoon")
    def read(n):
        G = s.amr(n, ["pressure", "alpha_v"])
        top = max(g["level"] for g in G)
        out = {"points": [], "pressure": [], "alpha_v": []}
        for g in G:
            if g["level"] != top:
                continue
            c = [0.5*(g[d][1:] + g[d][:-1]) for d in ("x", "y", "z")]
            Z, Y, X = np.meshgrid(c[2], c[1], c[0], indexing="ij")
            out["points"].append(np.stack([X.ravel(), Y.ravel(), Z.ravel()], 1))
            out["pressure"].append(g["pressure"].ravel())
            out["alpha_v"].append(g["alpha_v"].ravel())
        return {k: np.concatenate(v) for k, v in out.items()}, top
    return [(n, t) for n, t in enumerate(s.amr_times())], read

ap = argparse.ArgumentParser()
ap.add_argument("runs", nargs="*", default=["."])
ap.add_argument("--t0", type=float, default=None)
ap.add_argument("--t1", type=float, default=None)
ap.add_argument("--exp", default=None)
ap.add_argument("--plot", default=None)
args = ap.parse_args()

R = 0.5


def keyword(fn, key, n=1):
    for l in open(fn):
        s = l.split()
        if len(s) >= 2 + n and s[0] == key.split()[0] and s[1] == key.split()[1]:
            return [float(x) for x in s[2:2+n]]
    return None


results = []
for run in args.runs:
    rho = (keyword(f"{run}/ctrl.txt", "W 1") or [998.2])[0]
    W10 = keyword(f"{run}/ctrl.txt", "W 10")[0]
    b10 = keyword(f"{run}/control.txt", "B 10", 6)
    U = W10/((b10[3]-b10[2])*(b10[5]-b10[4]))
    q = 0.5*rho*U*U
    dx = keyword(f"{run}/control.txt", "B 127", 1)[0]
    # bins of the grid spacing, so that every bin holds band points on the cylinder
    bins = np.arange(0.0, 3.0 + 0.5*dx, dx)
    mid = 0.5*(bins[1:] + bins[:-1])
    cav = keyword(f"{run}/ctrl.txt", "W 60")
    pv = keyword(f"{run}/ctrl.txt", "W 61")[0] if cav else None

    amr = os.path.isdir(f"{run}/REEF3D_CFD.lagoon/amr/cfd_amr")
    if amr:
        st, read_amr = store_frames(run)
    else:
        st = steps(f"{run}/REEF3D_CFD_VTU")
    tmax = st[-1][1]
    t0 = args.t0 if args.t0 is not None else tmax - 5.0
    t1 = args.t1 if args.t1 is not None else tmax
    st = [(n, t) for n, t in st if t0 - 1e-6 <= t <= t1 + 1e-6]

    pr = np.loadtxt(glob.glob(f"{run}/REEF3D_CFD_PressureProbe/*-1.dat")[0], skiprows=6)
    w = (pr[:, 0] >= t0) & (pr[:, 0] <= t1)
    pref = pr[w, 1].mean()

    # band 1 (0.5-1.5 dx) and band 2 (1.5-2.5 dx): sums of p, distance, alpha and counts
    sp = np.zeros((2, len(mid))); sd = np.zeros((2, len(mid))); sa = np.zeros(len(mid)); cnt = np.zeros((2, len(mid)))
    dxb = dx
    for n, t in st:
        if amr:
            d, top = read_amr(n)
            dxb = dx/2**top      # bands at the spacing of the finest patches
        else:
            d = read_step(f"{run}/REEF3D_CFD_VTU", n)
        P = d["points"]
        x = P[:, 0]; r = np.hypot(P[:, 1], P[:, 2])
        nose = x < 0.0
        dist = np.where(nose, np.hypot(x, r) - R, r - R)
        s = np.where(nose, R*np.arctan2(r, -x), 0.25*math.pi + x)
        for b in (0, 1):
            m = (dist > (0.5+b)*dxb) & (dist < (1.5+b)*dxb) & (x < 7.0)
            k = np.digitize(s[m], bins) - 1
            ok = (k >= 0) & (k < len(mid))
            np.add.at(sp[b], k[ok], d["pressure"][m][ok])
            np.add.at(sd[b], k[ok], dist[m][ok])
            np.add.at(cnt[b], k[ok], 1)
            if b == 0:
                np.add.at(sa, k[ok], d["alpha_v"][m][ok] if "alpha_v" in d else 0.0)
    good = (cnt[0] > 0) & (cnt[1] > 0)
    Cp = np.full(len(mid), np.nan); al = np.full(len(mid), np.nan); Cp1 = np.full(len(mid), np.nan)
    p1 = sp[0, good]/cnt[0, good]; p2 = sp[1, good]/cnt[1, good]
    d1 = sd[0, good]/cnt[0, good]; d2 = sd[1, good]/cnt[1, good]
    Cp[good] = (p1 - (p2 - p1)*d1/(d2 - d1) - pref)/q
    Cp1[good] = (p1 - pref)/q
    al[good] = sa[good]/cnt[0, good]
    sigma = (pref - pv)/q if pv is not None else None

    # cavity: contiguous bins with alpha > 0.1 from the first one
    lcav = None
    c = np.where(al > 0.1)[0]
    if len(c):
        i = c[0]
        while i + 1 < len(mid) and al[i+1] > 0.1:
            i += 1
        lcav = (mid[c[0]], mid[i])
    imin = np.nanargmin(Cp)
    name = os.path.basename(os.path.abspath(run))
    print(f"{name}: window {t0:.2f}-{t1:.2f} s ({len(st)} frames{', AMR dx %.4g' % dxb if amr else ''}), U = {U:.3f}, p_ref = {pref:.2f} Pa, "
          + (f"sigma = {sigma:.3f}" if sigma is not None else "no cavitation")
          + f", Cp_min = {Cp[imin]:.3f} at s/D = {mid[imin]:.3f} (first band {np.nanmin(Cp1):.3f})"
          + (f", cavity s/D = {lcav[0]:.3f} - {lcav[1]:.3f}" if lcav else ""))
    results.append((name, sigma, mid, Cp, al))

print("\n  s/D  " + "".join(f"{n[:12]:>13}" for n, _, _, _, _ in results))
for sq in np.arange(0.1, 3.0, 0.1):
    print(f"{sq:6.2f} " + "".join(f"{np.interp(sq, m, Cp):13.3f}" for _, _, m, Cp, _ in results))

if args.plot:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    exp = np.loadtxt(args.exp, delimiter=",", skiprows=1) if args.exp else None
    fig, (a, b) = plt.subplots(2, 1, figsize=(8, 7), sharex=True, gridspec_kw={"height_ratios": [3, 1]})
    for (name, sg, mid, Cp, al), col in zip(results, plt.rcParams["axes.prop_cycle"].by_key()["color"]):
        lab = name if sg is None else f"{name}, sigma = {sg:.2f}"
        a.plot(mid, -Cp, color=col, lw=1.8, label=lab)
        if sg is not None:
            a.axhline(sg, color=col, lw=0.8, ls=":")
            b.plot(mid, al, color=col, lw=1.5)
        if exp is not None and len(exp):
            ref = sg if sg is not None else 99.0
            sel = np.abs(exp[:, 0] - ref) == np.abs(exp[:, 0] - ref).min()
            a.plot(exp[sel, 1], -exp[sel, 2], "o", color=col, mfc="none", ms=6)
    a.axvline(0.25*math.pi, color="k", lw=0.8, ls="--")
    a.set_ylabel("-Cp")
    a.legend(fontsize=9)
    a.grid(alpha=0.3)
    b.set_ylabel("alpha_v at the wall")
    b.set_xlabel("s/D from the nose tip (shoulder at 0.785)")
    b.set_xlim(0, 3)
    b.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(args.plot, dpi=130)
    print("written", os.path.abspath(args.plot))
