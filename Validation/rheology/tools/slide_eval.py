#!/usr/bin/env python3
# Architect: Hans Bihs
"""Evaluation of the two-layer landslide cases (Validation/rheology 16-19), SFLOW slide layer W 150.

  slide_eval.py mangeney <run slide> <run single layer or -> <out dir>   16: slide alone, Coulomb dam break vs Mangeney et al. (2000)
  slide_eval.py arrest <run slide> <run single layer or -> <out dir>     16: slide alone, viscoplastic arrest vs h^2 profile
  slide_eval.py rest <run> <out dir>                                     17: lake at rest over a slide at rest (2D)
  slide_eval.py subaerial <run> <out dir>                                17: subaerial slide into a lake (2D), volumes
  slide_eval.py dambreak <run> <out dir>                                 18: submarine dam break vs Ritter with g' = (1-r) g
  slide_eval.py slope35 <run> <out dir>                                  18: submarine granular slide on 35 deg, gauges

Reads the SFLOW free-surface VTP files (REEF3D_SFLOW_VTP_FSF, P 18 2), the slide log
(REEF3D_SFLOW_Log/REEF3D_SFLOW_slide.dat) and the wave gauges (P 51). Needs numpy and matplotlib.
"""
import sys, os, re, glob, math, json
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "turbulence", "tools"))
from r3read import read_vtk

BLUE, LBLUE, ORANGE, GREEN, GREY = "#2B4C7E", "#4F8FCF", "#E07B39", "#3A9A5B", "#6B6B6B"
KEYS = ("points", "velocity", "waterlevel", "eta", "wetdry", "slide_h", "slide_velocity", "slide_lambda")


def style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Liberation Sans", "Helvetica", "Arial"],
                         "font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.grid": True, "grid.color": "#DDDDDD", "grid.linewidth": 0.6, "axes.edgecolor": GREY,
                         "legend.frameon": False, "savefig.dpi": 200})
    return plt


def frames(run):
    files = sorted(glob.glob(os.path.join(run, "REEF3D_SFLOW_VTP_FSF", "REEF3D-SFLOW-FSF-*.vtp")))
    byn = {}
    for f in files:
        m = re.search(r"FSF-(\d+)-(\d+)\.vtp$", f)
        byn.setdefault(int(m.group(1)), []).append(f)
    out = []
    for n in sorted(byn):
        parts = [read_vtk(f) for f in byn[n]]
        head = open(byn[n][0], "rb").read(4000).decode("latin-1")
        t = float(re.search(r'Name="TimeValue"[^>]*>\s*([-0-9.eE+]+)', head).group(1))
        d = {k: np.concatenate([p[k] for p in parts]) for k in KEYS if k in parts[0]}
        out.append((t, d))
    return out


def at(F, t):
    return min(F, key=lambda f: abs(f[0] - t))


def line(d, key, y0=None, comp=None):
    """values of key along x at y closest to y0 (lowest y for 1D), duplicates removed"""
    P = d["points"]; x, y = P[:, 0], P[:, 1]
    ys = np.unique(np.round(y, 9))
    yb = ys[np.argmin(abs(ys - y0))] if y0 is not None else ys[0]
    sel = abs(y - yb) < 1e-9
    xs, idx = np.unique(np.round(x[sel], 9), return_index=True)
    ii = np.where(sel)[0][idx]
    v = d[key][ii]
    if comp is not None:
        v = v[:, comp]
    return xs, v


def front(x, h, thr):
    w = np.where(h > thr)[0]
    return x[w[-1]] if len(w) else np.nan


def slide_log(run):
    return np.loadtxt(os.path.join(run, "REEF3D_SFLOW_Log", "REEF3D_SFLOW_slide.dat"), comments="#")


def gauges(run):
    f = os.path.join(run, "REEF3D_SFLOW_WSF", "REEF3D-SFLOW-WSF-HG.dat")
    rows, xs = [], []
    for ln in open(f):
        s = ln.split()
        if not s:
            continue
        try:
            v = [float(a) for a in s]
        except ValueError:
            continue
        if len(v) == 3 and v[0] == int(v[0]) and not rows:
            xs.append(v[1]); continue
        rows.append(v)
    return np.array(xs), np.array([r for r in rows if len(r) == len(rows[-1])])

# ---------------------------------------------------------------------------------------------
def mangeney_exact(x, t, h0, gn, gx, tand):
    c0 = math.sqrt(gn*h0); m = gx - gn*tand
    xa = -c0*t + 0.5*m*t*t; xb = 2*c0*t + 0.5*m*t*t
    h = np.where(x <= xa, h0, 0.0); u = np.where(x <= xa, m*t, 0.0)
    mid = (x > xa) & (x < xb)
    h[mid] = (2*c0 - x[mid]/t + 0.5*m*t)**2/(9*gn)
    u[mid] = 2.0/3.0*(c0 + x[mid]/t + m*t)
    return h, u, xb


def mangeney(run, run1, out):
    plt = style()
    g = 9.81; th = math.radians(35); gx, gn = g*math.sin(th), g*math.cos(th); tand = math.tan(math.radians(25)); h0 = 1.0
    F = frames(run)
    F1 = frames(run1) if run1 not in (None, "-") and os.path.isdir(run1) else None
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
    rows = []
    cols = {1.0: LBLUE, 2.0: BLUE, 3.0: GREEN}
    for tt, col in cols.items():
        t, d = at(F, tt)
        x, h = line(d, "slide_h"); _, u = line(d, "slide_velocity", comp=0)
        he, ue, xb = mangeney_exact(x, t, h0, gn, gx, tand)
        c0 = math.sqrt(gn*h0); m = gx - gn*tand
        xr = -25.0 + c0*t + 0.5*m*t*t + 0.5
        sel = (x > xr) & (x < xb + 2)
        l1 = np.sum(abs(h - he)[sel])/np.sum(he[sel])
        wet = (h > 0.01) & (he > 0.01) & sel
        lu = np.sum(abs(u - ue)[wet])/np.sum(abs(ue[wet]))
        r = dict(t=round(t, 3), front_slide=round(float(front(x, h, 1e-3)), 3), front_exact=round(float(front(x, he, 1e-3)), 3),
                 L1_h=round(float(l1), 4), L1_u=round(float(lu), 4))
        if F1 is not None:
            t1, d1 = at(F1, tt)
            x1, h1 = line(d1, "waterlevel")
            hi = np.interp(x, x1, h1)
            r["front_single_layer"] = round(float(front(x1, h1, 1e-3)), 3)
            r["max_diff_to_single_layer"] = round(float(np.max(abs(h - hi))), 5)
            ax[0].plot(x1, h1, color=ORANGE, lw=1.0, ls=":")
        rows.append(r)
        ax[0].plot(x, he, color=GREY, lw=1.0, ls="--")
        ax[0].plot(x, h, color=col, lw=1.6, label=f"slide layer, t = {tt:.0f} s")
        ax[1].plot(x[he > 1e-3], ue[he > 1e-3], color=GREY, lw=1.0, ls="--")
        ax[1].plot(x[h > 1e-3], u[h > 1e-3], color=col, lw=1.6)
    ax[0].plot([], [], color=GREY, ls="--", lw=1.0, label="Mangeney et al. (2000)")
    if F1 is not None:
        ax[0].plot([], [], color=ORANGE, ls=":", lw=1.0, label="single layer W 120 (case 01)")
    for a in ax:
        a.set_xlim(-12, 30); a.set_xlabel("x along slope (m)")
    ax[0].set_ylabel("h (m)"); ax[0].legend(fontsize=8); ax[1].set_ylabel("u (m/s)")
    fig.tight_layout(); fig.savefig(os.path.join(out, "16_slide_mangeney.png"))
    L = slide_log(run)
    return dict(profiles=rows, max_rel_volume_error=float(np.max(abs(L[:, 4]))))


def arrest(run, run1, out):
    plt = style()
    g, rho, ty = 9.81, 1000.0, 89.0
    V = 0.5*0.3
    a = math.sqrt(2*ty/(rho*g)); xf_ex = (3*V/(2*a))**(2.0/3.0)
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
    res = {}
    for name, rr, key, col in (("slide layer W 150", run, "slide_h", BLUE), ("single layer W 120 (case 02)", run1, "waterlevel", ORANGE)):
        if rr in (None, "-") or not os.path.isdir(rr):
            continue
        F = frames(rr)
        ts, um = [], []
        for t, d in F:
            x, h = line(d, key)
            _, u = line(d, "slide_velocity" if key == "slide_h" else "velocity", comp=0)
            ts.append(t); um.append(np.max(abs(u[h > 1e-3])) if np.any(h > 1e-3) else 0.0)
        x, h = line(F[-1][1], key)
        res[name] = dict(t_end=round(F[-1][0], 2), front=round(float(front(x, h, 1e-3)), 4), h0=round(float(h[0]), 4),
                         umax_end=float(um[-1]), volume=round(float(np.trapezoid(h, x)), 5))
        ax[0].plot(x, h, color=col, lw=1.6 if key == "slide_h" else 1.0, ls="-" if key == "slide_h" else ":", label=name)
        ax[1].semilogy(ts, np.maximum(um, 1e-12), color=col, lw=1.6, label=name)
    xx = np.linspace(0, xf_ex, 400)
    ax[0].plot(xx, a*np.sqrt(xf_ex - xx), color=GREY, ls="--", lw=1.2, label=r"$h^2 = 2\tau_y(x_f-x)/\rho g$")
    ax[0].set_xlim(0, 2); ax[0].set_xlabel("x (m)"); ax[0].set_ylabel("h (m)"); ax[0].legend(fontsize=8)
    ax[1].set_xlabel("t (s)"); ax[1].set_ylabel("max |u| (m/s)"); ax[1].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(out, "16_slide_arrest.png"))
    res["exact_fully_yielded"] = dict(front=round(xf_ex, 4), h0=round(a*math.sqrt(xf_ex), 4))
    return res


def rest(run, out):
    plt = style()
    F = frames(run)
    L = slide_log(run)
    um_w, um_s, de = [], [], []
    eta0 = F[0][1]["eta"]; wet0 = F[0][1]["waterlevel"] > 1e-3
    for t, d in F:
        w = d["waterlevel"] > 1e-3
        um_w.append(np.max(np.linalg.norm(d["velocity"][w, :2], axis=1)) if np.any(w) else 0.0)
        s = d["slide_h"] > 1e-3
        um_s.append(np.max(np.linalg.norm(d["slide_velocity"][s, :2], axis=1)) if np.any(s) else 0.0)
        de.append(np.max(abs(d["eta"][w & wet0] - eta0[w & wet0])))
    vw = L[:, 5]
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.4))
    ts = [f[0] for f in F]
    ax[0].semilogy(ts, np.maximum(um_w, 1e-16), color=BLUE, lw=1.6, label="water max |u|")
    ax[0].semilogy(ts, np.maximum(um_s, 1e-16), color=ORANGE, lw=1.6, label="slide max |u|")
    ax[0].semilogy(ts, np.maximum(de, 1e-16), color=GREEN, lw=1.6, label=r"max |$\eta - \eta_0$|")
    ax[0].set_xlabel("t (s)"); ax[0].legend(fontsize=8)
    t, d = F[-1]
    x, zb = line(d, "slide_h", y0=1.0)
    _, eta = line(d, "eta", y0=1.0)
    ax[1].plot(x, eta + 1.0, color=BLUE, lw=1.4, label="water surface")
    ax[1].plot(x, zb, color=ORANGE, lw=1.4, label="slide depth")
    ax[1].set_xlabel("x (m) at y = 1 m"); ax[1].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(out, "17_rest.png"))
    return dict(t_end=round(t, 2), water_umax=float(max(um_w)), slide_umax=float(max(um_s)), eta_dev_max=float(max(de)),
                slide_rel_volume_error=float(np.max(abs(L[:, 4]))), water_rel_volume_error=float(np.max(abs(vw[1:] - vw[1]))/vw[1]))


def subaerial(run, out):
    plt = style()
    F = frames(run)
    L = slide_log(run)
    fig, ax = plt.subplots(2, 2, figsize=(10, 6.0))
    a = ax[0, 0]
    for tt, col in ((0.0, GREY), (1.0, LBLUE), (2.0, BLUE), (4.0, GREEN), (8.0, ORANGE)):
        t, d = at(F, tt)
        x, hs = line(d, "slide_h", y0=2.0)
        _, wl = line(d, "waterlevel", y0=2.0)
        zb = np.where(x > 12, 2.0*(x - 12)/8.0, 0.0)
        a.plot(x, zb + hs, color=col, lw=1.4, label=f"t = {t:.1f} s")
        a.plot(x, np.where(wl > 1e-3, zb + hs + wl, np.nan), color=col, lw=1.0, ls="--")
    a.plot(x, zb, color="k", lw=0.8)
    a.set_xlim(8, 20); a.set_ylim(-0.1, 1.6); a.set_xlabel("x (m) at y = 2 m"); a.set_ylabel("z (m)"); a.legend(fontsize=7)
    a = ax[0, 1]
    a.plot(L[:, 0], L[:, 4], color=ORANGE, lw=1.4, label="slide volume rel. error")
    a.plot(L[1:, 0], (L[1:, 5] - L[1, 5])/L[1, 5], color=BLUE, lw=1.4, label="water volume rel. error")
    a.set_xlabel("t (s)"); a.legend(fontsize=8)
    a = ax[1, 0]
    a.plot(L[:, 0], L[:, 7], color=ORANGE, lw=1.4, label="slide max |u|")
    a.plot(L[:, 0], L[:, 8], color=BLUE, lw=1.4, label="slide x_min (front)")
    a.plot(L[:, 0], L[:, 9], color=LBLUE, lw=1.4, label="slide x_max (rear)")
    a.set_xlabel("t (s)"); a.legend(fontsize=8)
    a = ax[1, 1]
    em = []
    for t, d in F:
        w = d["waterlevel"] > 1e-3
        em.append((t, np.max(d["eta"][w]), np.min(d["eta"][w])))
    em = np.array(em)
    a.plot(em[:, 0], em[:, 1], color=BLUE, lw=1.4, label=r"max $\eta$")
    a.plot(em[:, 0], em[:, 2], color=LBLUE, lw=1.4, label=r"min $\eta$")
    a.set_xlabel("t (s)"); a.set_ylabel(r"$\eta$ (m)"); a.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(out, "17_subaerial.png"))
    return dict(slide_rel_volume_error=float(np.max(abs(L[:, 4]))), water_rel_volume_error=float(np.max(abs(L[1:, 5] - L[1, 5]))/L[1, 5]),
                eta_max=float(np.max(em[:, 1])), eta_min=float(np.min(em[:, 2])), slide_umax=float(np.max(L[:, 7])),
                substeps_max=int(np.max(L[:, 2])))


def dambreak(run, out):
    plt = style()
    g, rhos, rhow, h0, x0 = 9.81, 2000.0, 998.2, 1.0, 10.0
    gp = g*(1 - rhow/rhos); c0 = math.sqrt(gp*h0)
    F = frames(run)
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
    rows = []
    for tt, col in ((1.0, LBLUE), (2.0, BLUE), (3.0, GREEN), (4.0, ORANGE)):
        t, d = at(F, tt)
        x, h = line(d, "slide_h"); _, u = line(d, "slide_velocity", comp=0); _, eta = line(d, "eta")
        xi = (x - x0)/t
        he = np.where(xi < -c0, h0, np.where(xi > 2*c0, 0.0, (2*c0 - xi)**2/(9*gp)))
        sel = x > x0 - c0*t + 0.5
        l1 = np.sum(abs(h - he)[sel])/np.sum(he[sel])
        rows.append(dict(t=round(t, 3), front=round(float(front(x, h, 1e-3)), 3), front_ritter=round(x0 + 2*c0*t, 3),
                         L1_h=round(float(l1), 4), eta_range=[round(float(eta.min()), 5), round(float(eta.max()), 5)]))
        ax[0].plot(x, he, color=GREY, lw=1.0, ls="--")
        ax[0].plot(x, h, color=col, lw=1.6, label=f"t = {tt:.0f} s")
        ax[1].plot(x, eta, color=col, lw=1.4)
    ax[0].plot([], [], color=GREY, ls="--", label=r"Ritter, $g' = (1-r)g$")
    ax[0].set_xlim(0, 35); ax[0].set_xlabel("x (m)"); ax[0].set_ylabel(r"$h_s$ (m)"); ax[0].legend(fontsize=8)
    ax[1].set_xlim(0, 60); ax[1].set_xlabel("x (m)"); ax[1].set_ylabel(r"$\eta$ (m)")
    fig.tight_layout(); fig.savefig(os.path.join(out, "18_dambreak.png"))
    L = slide_log(run)
    return dict(g_reduced=round(gp, 4), profiles=rows, slide_rel_volume_error=float(np.max(abs(L[:, 4]))),
                water_rel_volume_error=float(np.max(abs(L[1:, 5] - L[1, 5]))/L[1, 5]))


def slope35(run, out, run_b=None, label_a="hydrostatic", label_b=None):
    plt = style()
    F = frames(run)
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.8))
    a = ax[0]
    for tt, col in ((0.0, GREY), (0.3, LBLUE), (0.6, BLUE), (1.0, GREEN), (2.0, ORANGE)):
        t, d = at(F, tt)
        x, hs = line(d, "slide_h")
        zb = np.where(x < 0.6, 0.42 - 0.7*x, 0.0)
        a.plot(x, zb + hs, color=col, lw=1.4, label=f"t = {t:.1f} s")
        _, wl = line(d, "waterlevel"); _, eta = line(d, "eta")
        a.plot(x, np.where(wl > 1e-4, 0.33 + eta, np.nan), color=col, lw=0.9, ls="--")
    a.plot(x, zb, color="k", lw=0.8); a.set_xlim(0, 2.0); a.set_ylim(-0.01, 0.40)
    a.set_xlabel("x (m)"); a.set_ylabel("z (m)"); a.legend(fontsize=7)
    xs, G = gauges(run)
    for k in range(G.shape[1] - 1):
        ax[1].plot(G[:, 0], 100*G[:, k + 1] + 3*k, color=[BLUE, LBLUE, GREEN, ORANGE][k % 4], lw=1.3, label=f"x = {xs[k]:.1f} m" if len(xs) > k else None)
    ax[1].set_xlabel("t (s)"); ax[1].set_ylabel(r"$\eta$ (cm), offset 3 cm per gauge"); ax[1].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(out, "18_slope35.png"))
    L = slide_log(run)
    return dict(slide_rel_volume_error=float(np.max(abs(L[:, 4]))), water_rel_volume_error=float(np.max(abs(L[1:, 5] - L[1, 5]))/L[1, 5]),
                gauge_max_cm=[round(100*float(np.max(G[:, k + 1])), 3) for k in range(G.shape[1] - 1)],
                gauge_min_cm=[round(100*float(np.min(G[:, k + 1])), 3) for k in range(G.shape[1] - 1)],
                slide_umax=float(np.max(L[:, 7])), runout=float(np.max(L[:, 9])))


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd in ("mangeney", "arrest"):
        r = {"mangeney": mangeney, "arrest": arrest}[cmd](sys.argv[2], sys.argv[3], sys.argv[4])
    else:
        r = {"rest": rest, "subaerial": subaerial, "dambreak": dambreak, "slope35": slope35}[cmd](sys.argv[2], sys.argv[3])
    print(json.dumps(r, indent=1))
