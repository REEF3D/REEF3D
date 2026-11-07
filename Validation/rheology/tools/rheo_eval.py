#!/usr/bin/env python3
# Architect: Hans Bihs
"""Evaluation of the rheology validation cases (Validation/rheology), SFLOW rheology core W 120.

  rheo_eval.py mangeney <run dir> <out dir>      01: Coulomb dam break on an incline vs Mangeney et al. (2000)
  rheo_eval.py arrest <run core> <run legacy> <out dir>   02: viscoplastic arrest vs h^2 = 2 tau_y (x_f - x)/(rho g)
  rheo_eval.py cochard <run dir> <out dir>       03: Carbopol 0.30 %, 12 deg vs Cochard & Ancey (2009)
  rheo_eval.py pore <run dir> <out dir>          04: uniform layer with dissipating pore pressure vs analytic u(t)
  rheo_eval.py uniform <runs prefix> <ode binary> <out dir>   05: uniform layers, all basal laws, vs the ODE u(t)

Reads the SFLOW free-surface VTP files (REEF3D_SFLOW_VTP_FSF); needs numpy and matplotlib.
"""
import sys, os, re, glob, math, json
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "turbulence", "tools"))
from r3read import read_vtk

BLUE, LBLUE, ORANGE, GREEN, GREY = "#2B4C7E", "#4F8FCF", "#E07B39", "#3A9A5B", "#6B6B6B"

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
    """list of (time, dict) of the merged pieces of every printout, sorted by time"""
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
        d = {}
        for k in ("points", "velocity", "waterlevel", "eta", "wetdry"):
            if k in parts[0]:
                d[k] = np.concatenate([p[k] for p in parts])
        out.append((t, d))
    return out

def line(d, y0=None, tol=None):
    """points near y = y0 (all for 1D), sorted in x, duplicates (piece borders) removed"""
    P = d["points"]; x, y = P[:, 0], P[:, 1]
    sel = np.ones(len(x), bool)
    if y0 is not None:
        ys = np.unique(np.round(y, 9))
        yb = ys[np.argmin(abs(ys - y0))]
        sel = abs(y - yb) < 1e-9
    else:
        yb = np.min(y)
        sel = abs(y - yb) < 1e-9
    xs, idx = np.unique(np.round(x[sel], 9), return_index=True)
    ii = np.where(sel)[0][idx]
    return xs, d["waterlevel"][ii], d["velocity"][ii, 0]

def front(x, h, thr):
    w = np.where(h > thr)[0]
    return x[w[-1]] if len(w) else np.nan

# ---------------------------------------------------------------------------------------------
def mangeney_exact(x, t, h0, gn, gx, tand):
    c0 = math.sqrt(gn*h0); m = gx - gn*tand
    xa = -c0*t + 0.5*m*t*t; xb = 2*c0*t + 0.5*m*t*t
    h = np.where(x <= xa, h0, 0.0); u = np.where(x <= xa, m*t, 0.0)
    mid = (x > xa) & (x < xb)
    h[mid] = (2*c0 - x[mid]/t + 0.5*m*t)**2/(9*gn)
    u[mid] = 2.0/3.0*(c0 + x[mid]/t + m*t)
    return h, u, xb

def mangeney(run, out):
    plt = style()
    g = 9.81; th = math.radians(35); gx, gn = g*math.sin(th), g*math.cos(th); tand = math.tan(math.radians(25)); h0 = 1.0
    F = frames(run)
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
    rows = []
    cols = {1.0: LBLUE, 2.0: BLUE, 3.0: GREEN}
    for t, d in F:
        tt = round(t*2)/2
        if tt not in cols or abs(t - tt) > 0.06:
            continue
        x, h, u = line(d)
        he, ue, xb = mangeney_exact(x, t, h0, gn, gx, tand)
        xf = front(x, h, 1e-3)
        xfe = front(x, he, 1e-3)
        c0 = math.sqrt(gn*h0); m = gx - gn*tand
        xr = -25.0 + c0*t + 0.5*m*t*t + 0.5        # head of the rarefaction from the upstream wall
        sel = (x > xr) & (x < xb + 2)
        l1 = np.sum(abs(h - he)[sel])/np.sum(he[sel])
        wet = (h > 0.01) & (he > 0.01) & sel
        lu = np.sum(abs(u - ue)[wet])/np.sum(abs(ue[wet]))
        rows.append(dict(t=round(t, 3), front_model=round(float(xf), 3), front_exact=round(float(xfe), 3), front_exact_h0=round(float(xb), 3),
                         L1_h=round(float(l1), 4), L1_u=round(float(lu), 4)))
        ax[0].plot(x, he, color=GREY, lw=1.0, ls="--")
        ax[0].plot(x, h, color=cols[tt], lw=1.6, label=f"t = {tt:.0f} s")
        ax[1].plot(x[he > 1e-3], ue[he > 1e-3], color=GREY, lw=1.0, ls="--")
        ax[1].plot(x[h > 1e-3], u[h > 1e-3], color=cols[tt], lw=1.6)
    ax[0].plot([], [], color=GREY, ls="--", lw=1.0, label="Mangeney et al. (2000)")
    ax[0].set_xlim(-12, 30); ax[1].set_xlim(-12, 30)
    ax[0].set_xlabel("x along slope (m)"); ax[0].set_ylabel("h (m)"); ax[0].legend()
    ax[1].set_xlabel("x along slope (m)"); ax[1].set_ylabel("u (m/s)")
    fig.tight_layout(); fig.savefig(os.path.join(out, "01_mangeney_profiles.png")); 
    return rows

# ---------------------------------------------------------------------------------------------
def arrest(run_core, run_legacy, out):
    plt = style()
    g, rho, ty = 9.81, 1000.0, 89.0
    V = 0.5*0.3
    a = math.sqrt(2*ty/(rho*g)); xf_ex = (3*V/(2*a))**(2.0/3.0)
    res = {}
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
    for name, run, col in (("rheology core W 120 5", run_core, BLUE), ("legacy W 90 1", run_legacy, ORANGE)):
        if run is None or not os.path.isdir(run):
            continue
        F = frames(run)
        ts, xf, um = [], [], []
        for t, d in F:
            x, h, u = line(d)
            ts.append(t); xf.append(front(x, h, 1e-3)); um.append(np.max(abs(u[h > 1e-3])) if np.any(h > 1e-3) else 0.0)
        x, h, u = line(F[-1][1])
        res[name] = dict(t_end=round(F[-1][0], 2), front=round(float(xf[-1]), 4), umax_end=float(um[-1]),
                         volume=round(float(np.trapezoid(h, x)), 5))
        ax[0].plot(x, h, color=col, lw=1.6, label=name)
        ax[1].semilogy(ts, np.maximum(um, 1e-12), color=col, lw=1.6, label=name)
    xx = np.linspace(0, xf_ex, 400)
    ax[0].plot(xx, a*np.sqrt(xf_ex - xx), color=GREY, ls="--", lw=1.2, label=r"$h^2 = 2\tau_y(x_f-x)/\rho g$")
    ax[0].axhline(0.3, xmax=0.5/3.0, color=LBLUE, lw=1.0, ls=":")
    ax[0].set_xlim(0, 3); ax[0].set_xlabel("x (m)"); ax[0].set_ylabel("h (m)"); ax[0].legend()
    ax[1].set_xlabel("t (s)"); ax[1].set_ylabel("max |u| (m/s)"); ax[1].legend()
    fig.tight_layout(); fig.savefig(os.path.join(out, "02_arrest.png"))
    res["exact"] = dict(front=round(xf_ex, 4), h_wall=round(a*math.sqrt(xf_ex), 4))
    return res

# ---------------------------------------------------------------------------------------------
def cochard(run, out):
    plt = style()
    F = frames(run)
    ts, xf = [], []
    for t, d in F:
        x, h, u = line(d, y0=0.0)
        f = front(x, h, 1e-3)
        ts.append(t); xf.append(f)
    ts = np.array(ts); xf = np.array(xf)
    # fit of the late regime
    sel = ts >= 5.0
    p = np.polyfit(np.log(ts[sel]), np.log(xf[sel]), 1) if np.sum(sel) > 3 else [np.nan]
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
    ax[0].plot(ts, xf, color=BLUE, lw=1.6, label="SFLOW W 120 5")
    ax[0].axhspan(1.20, 1.45, color=ORANGE, alpha=0.2, lw=0, label="C&A: end of inertial phase (1.20-1.45 m)")
    ax[0].plot([0, 1.0], [0, 0.7], color=ORANGE, lw=1.2, ls="--", label="C&A: front speed ~0.7 m/s, t < 1 s")
    ax[0].set_xlabel("t (s)"); ax[0].set_ylabel("front position x_f (m)"); ax[0].legend(loc="lower right")
    ax[1].loglog(ts[ts > 0], xf[ts > 0], color=BLUE, lw=1.6, label=f"SFLOW, late slope {p[0]:.3f}")
    tt = np.array([5, ts.max()]); ax[1].loglog(tt, xf[sel][0]*(tt/tt[0])**0.275, color=ORANGE, ls="--", lw=1.2, label=r"C&A: $x_f \propto t^{0.275}$")
    ax[1].set_xlabel("t (s)"); ax[1].set_ylabel("x_f (m)"); ax[1].legend()
    fig.tight_layout(); fig.savefig(os.path.join(out, "03_cochard_front.png"))
    # map of the deposit
    t, d = F[-1]
    P = d["points"]
    fig2, a2 = plt.subplots(figsize=(6, 3.6))
    hh = d["waterlevel"]
    xu = np.unique(np.round(P[:, 0], 9)); yu = np.unique(np.round(P[:, 1], 9))
    Hm = np.full((len(yu), len(xu)), np.nan)
    ix = np.searchsorted(xu, np.round(P[:, 0], 9)); iy = np.searchsorted(yu, np.round(P[:, 1], 9))
    Hm[iy, ix] = np.where(hh > 1e-3, hh, np.nan)
    m = a2.pcolormesh(xu, yu, Hm, cmap="Blues", shading="nearest", vmin=0)
    a2.plot([-0.52, 0, 0, -0.52], [0.15, 0.15, -0.15, -0.15], color=GREY, lw=1.0)
    a2.set_aspect("equal"); a2.set_xlabel("x (m)"); a2.set_ylabel("y (m)"); a2.set_title(f"h at t = {t:.1f} s", fontsize=10)
    fig2.colorbar(m, ax=a2, label="h (m)"); fig2.tight_layout(); fig2.savefig(os.path.join(out, "03_cochard_map.png"))
    k1 = np.argmin(abs(ts - 1.0)); k2 = np.argmin(abs(ts - 2.0))
    return dict(front_1s=round(float(xf[k1]), 3), front_2s=round(float(xf[k2]), 3), front_end=round(float(xf[-1]), 3),
                t_end=round(float(ts[-1]), 2), mean_speed_0_1s=round(float(xf[k1]/ts[k1]), 3), late_slope=round(float(p[0]), 3))

# ---------------------------------------------------------------------------------------------
def pore(run, out):
    plt = style()
    g = 9.81; th = math.radians(20); gx, gn = g*math.sin(th), g*math.cos(th); tand = math.tan(math.radians(25))
    l0, lh, h = 0.5, 0.0, 1.0; D = 4*h*h/(math.pi**2*5.0); T = 4*h*h/(math.pi**2*D)
    A = gx - gn*tand*(1 - lh); B = gn*tand*(l0 - lh)
    def uex(t):
        u = A*t + B*T*(1 - np.exp(-t/T))
        # after the first zero crossing (t > 0) the layer stays at rest
        tt = np.linspace(1e-6, 60, 600001); uu = A*tt + B*T*(1 - np.exp(-tt/T)); ts = tt[np.argmax(uu < 0)]
        return np.where(t < ts, u, 0.0), ts
    F = frames(run)
    ts, um = [], []
    for t, d in F:
        x, hh, u = line(d)
        k = np.argmin(abs(x - 150.0)); ts.append(t); um.append(u[k])
    ts = np.array(ts); um = np.array(um)
    ue, tstop = uex(ts)
    tf = np.linspace(0, ts.max(), 600); uf, _ = uex(tf)
    fig, ax = plt.subplots(figsize=(5.5, 3.6))
    ax.plot(tf, uf, color=GREY, ls="--", lw=1.2, label="analytic")
    ax.plot(ts, um, "o", color=BLUE, ms=3.5, label="SFLOW W 120 1, W 126 2")
    ax.axvline(tstop, color=ORANGE, lw=1.0, ls=":")
    ax.set_xlabel("t (s)"); ax.set_ylabel("u at x = 150 m (m/s)"); ax.legend()
    fig.tight_layout(); fig.savefig(os.path.join(out, "04_pore_pressure.png"))
    move = um > 1e-6
    t_stop_model = ts[np.argmax((ts > 1.0) & ~move)] if np.any((ts > 1.0) & ~move) else np.nan
    return dict(max_abs_err=round(float(np.max(abs(um - ue))), 4), u_peak_model=round(float(um.max()), 4),
                u_peak_exact=round(float(uf.max()), 4), t_stop_exact=round(float(tstop), 3),
                t_stop_model_first_output_at_rest=float(t_stop_model), u_end=float(um[-1]))

# ---------------------------------------------------------------------------------------------
UNIFORM = [  # name, law, theta, h, rho, tend, ODE parameters, label, centre x
    ("voellmy", 2, 30, 1.0, 2000.0, 8.0, [21.801, 200], "Voellmy", 150.0),
    ("pouliquen", 3, 25, 0.005, 1500.0, 2.0, [21.0, 30.7, 22.2, 0.136, 0.00065], "Pouliquen-Forterre", 2.0),
    ("muI", 4, 25, 0.005, 1500.0, 2.0, [20.9, 32.76, 0.279, 0.0005, 0.6], "mu(I)", 2.0),
    ("hb", 5, 12, 0.1, 1000.0, 3.0, [89.0, 47.68, 0.415, 0.0], "Herschel-Bulkley", 5.0),
    ("hbc", 6, 12, 0.1, 1000.0, 3.0, [20.0, 47.68, 0.415, 5.711], "Herschel-Bulkley-Coulomb", 5.0),
]

def uniform(prefix, ode, out):
    import subprocess
    plt = style()
    fig, axs = plt.subplots(1, 5, figsize=(15, 3.2))
    rows = {}
    for ax, (n, law, th, h, rho, tend, par, lab, xc) in zip(axs, UNIFORM):
        txt = subprocess.run([ode, str(law), str(th), str(h), str(rho), str(tend)] + [str(v) for v in par],
                             capture_output=True, text=True).stdout.split()
        te = np.array(txt[0::2], float); ue = np.array(txt[1::2], float)
        ts, um = [], []
        for t, d in frames(prefix + n):
            x, hh, u = line(d); k = np.argmin(abs(x - xc)); ts.append(t); um.append(u[k])
        ts = np.array(ts); um = np.array(um); ui = np.interp(ts, te, ue)
        rows[n] = dict(u_end_model=round(float(um[-1]), 5), u_end_ode=round(float(ui[-1]), 5),
                       max_rel_err=round(float(np.max(abs(um - ui))/np.max(ue)), 4))
        ax.plot(te, ue, color=GREY, ls="--", lw=1.2, label="ODE (core)")
        ax.plot(ts, um, "o", color=BLUE, ms=2.5, label="SFLOW")
        ax.set_title(lab, fontsize=10); ax.set_xlabel("t (s)")
    axs[0].set_ylabel("u (m/s)"); axs[0].legend()
    fig.tight_layout(); fig.savefig(os.path.join(out, "05_uniform_layers.png"))
    return rows

if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "mangeney":
        r = mangeney(sys.argv[2], sys.argv[3])
    elif cmd == "arrest":
        r = arrest(sys.argv[2], sys.argv[3], sys.argv[4])
    elif cmd == "cochard":
        r = cochard(sys.argv[2], sys.argv[3])
    elif cmd == "pore":
        r = pore(sys.argv[2], sys.argv[3])
    elif cmd == "uniform":
        r = uniform(sys.argv[2], sys.argv[3], sys.argv[4])
    print(json.dumps(r, indent=1))
