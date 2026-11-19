#!/usr/bin/env python3
# Architect: Hans Bihs
"""Evaluation of the rheology validation cases (Validation/rheology), SFLOW rheology core W 120.

  rheo_eval.py mangeney <run dir> <out dir>      01: Coulomb dam break on an incline vs Mangeney et al. (2000)
  rheo_eval.py arrest <run core> <run legacy> <out dir>   02: viscoplastic arrest vs h^2 = 2 tau_y (x_f - x)/(rho g)
  rheo_eval.py cochard <run dir> <out dir>       03: Carbopol 0.30 %, 12 deg vs Cochard & Ancey (2009)
  rheo_eval.py pore <run dir> <out dir>          04: uniform layer with dissipating pore pressure vs analytic u(t)
  rheo_eval.py uniform <runs prefix> <ode binary> <out dir>   05: uniform layers, all basal laws, vs the ODE u(t)
  rheo_eval.py films <runs prefix> <out dir>     06: CFD steady films on an incline (HB, HBC, mu(I)) vs the exact profiles
  rheo_eval.py collapse <runs prefix> <out dir>  07: CFD mu(I) column collapse a = 0.5, 2 vs Lube et al. (2005)
  rheo_eval.py cfdarrest <run cfd> <run sflow> <out dir>   08: CFD viscoplastic arrest vs h^2 profile and SFLOW (02)
  rheo_eval.py terrainslide <run> <case> <out> <name>     26: real landslide on DEMs (slide_terrain.py) vs the DEM deposit
  rheo_eval.py terrainseries <case> <out> <name>:<value>:<run> ...   26: a parameter series of terrainslide runs
  rheo_eval.py tsunamirunup <run> <case> <runup.csv> <out> <name> [radius]   27: landslide tsunami runup vs a survey
  rheo_eval.py tsunamiseries <case> <runup.csv> <out> <name>:<run> ...       27: a series of tsunamirunup runs

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

# ---------------------------------------------------------------------------------------------
# CFD
def cfd_frames(run, last_only=False):
    files = sorted(glob.glob(os.path.join(run, "REEF3D_CFD_VTU", "*.vtu")))
    byn = {}
    for f in files:
        m = re.search(r"-(\d+)-(\d+)\.vtu$", f)
        if m:
            byn.setdefault(int(m.group(1)), []).append(f)
    keys = sorted(byn)
    if last_only:
        keys = keys[-1:]
    out = []
    for n in keys:
        parts = [read_vtk(f) for f in byn[n]]
        head = open(byn[n][0], "rb").read(4000).decode("latin-1")
        mt = re.search(r'Name="TimeValue"[^>]*>\s*([-0-9.eE+]+)', head)
        t = float(mt.group(1)) if mt else float(n)
        d = {k: np.concatenate([p[k] for p in parts]) for k in ("points", "velocity", "phi") if k in parts[0]}
        out.append((t, d))
    return out

FILM_DELTA = {"hb": 0.0075, "hbc": 0.0, "muI": 0.0}   # hb: lithostatic pressure from the level set above the lid (no effect), hbc, muI: solved pressure (W 131 2)

def film_exact(law, z, h, th, rho, par):
    g = 9.81
    dl = FILM_DELTA[law]
    if law in ("hb", "hbc"):
        tauy, K, n = par[0], par[1], par[2]
        G = rho*g*math.sin(th) if law == "hb" else rho*g*(math.sin(th) - math.tan(math.radians(par[3]))*math.cos(th))
        if law == "hbc":
            tauy = tauy + math.tan(math.radians(par[3]))*rho*g*math.cos(th)*dl
        h0 = h - tauy/G
        a = (G/K)**(1.0/n)
        zz = np.minimum(z, h0)
        return n/(n + 1.0)*a*(h0**((n + 1.0)/n) - (h0 - zz)**((n + 1.0)/n))
    # mu(I) Bagnold profile: mu(I) = tan(theta), I constant
    mus, mu2, I0, d, phis = math.tan(math.radians(par[0])), math.tan(math.radians(par[1])), par[2], par[3], par[4]
    rhos = rho/phis
    # mu(I(z)) = tau/p' = tan(theta) (h - z)/(h + delta - z), gamma = I sqrt(p'/rho_s)/d, u = int gamma dz
    zf = np.linspace(0.0, h, 4001)
    r = math.tan(th)*(h - zf)/(h + dl - zf)
    I = np.where(r > mus, I0*(r - mus)/np.maximum(mu2 - r, 1e-12), 0.0)
    gam = I*np.sqrt(rho*g*math.cos(th)*(h + dl - zf)/rhos)/d
    uf = np.concatenate([[0.0], np.cumsum(0.5*(gam[1:] + gam[:-1])*np.diff(zf))])
    return np.interp(z, zf, uf)

FILMS = [("hb", 12, 0.1, 1000.0, [89.0, 47.68, 0.415], "Herschel-Bulkley"),
         ("hbc", 12, 0.1, 1000.0, [20.0, 47.68, 0.415, 5.711], "Herschel-Bulkley-Coulomb"),
         ("muI", 25, 0.02, 1500.0, [20.9, 32.76, 0.279, 0.0005, 0.6], "mu(I) (Bagnold)")]

def d_pressure_line(run, sel, o):
    t, d = cfd_frames(run, last_only=True)[0]
    files = sorted(glob.glob(os.path.join(run, "REEF3D_CFD_VTU", "*.vtu")))
    p = np.concatenate([read_vtk(f)["pressure"] for f in files if re.search(r"-%08d-" % int(re.search(r"-(\d+)-\d+\.vtu$", files[-1]).group(1)), f)])
    return p[sel][o]

def film_ref_p(law, z, p, h, th, rho, par):
    """exact steady profile for the pressure p(z) of the run: tau = rho g sin(theta) (h - z), p' = p(z)"""
    g = 9.81
    zf = np.linspace(0.0, h, 4001); pf = np.maximum(np.interp(zf, z, p), 0.0)
    tau = rho*g*math.sin(th)*(h - zf)
    if law == "hbc":
        tauy = par[0] + math.tan(math.radians(par[3]))*pf
        gam = np.where(tau > tauy, (np.maximum(tau - tauy, 0.0)/par[1])**(1.0/par[2]), 0.0)
    else:
        mus, mu2, I0, d, phis = math.tan(math.radians(par[0])), math.tan(math.radians(par[1])), par[2], par[3], par[4]
        r = np.where(pf > 0, tau/np.maximum(pf, 1e-12), mu2)
        r = np.minimum(r, mu2 - 1e-3)
        I = np.where(r > mus, I0*(r - mus)/(mu2 - r), 0.0)
        gam = I*np.sqrt(pf/(rho/phis))/d
    uf = np.concatenate([[0.0], np.cumsum(0.5*(gam[1:] + gam[:-1])*np.diff(zf))])
    return np.interp(z, zf, uf)

def films(prefix, out):
    plt = style()
    fig, axs = plt.subplots(1, 3, figsize=(11, 3.6))
    rows = {}
    for ax, (n, th, h, rho, par, lab) in zip(axs, FILMS):
        F = cfd_frames(prefix + n)
        t, d = F[-1]
        P = d["points"]; x = P[:, 0]; z = P[:, 2]
        xs = np.unique(np.round(x, 9)); xm = xs[len(xs)//2]
        sel = abs(x - xm) < 1e-9
        zz = z[sel]; uu = d["velocity"][sel, 0]; ph = d["phi"][sel]
        o = np.argsort(zz); zz, uu, ph = zz[o], uu[o], ph[o]
        pp = d_pressure_line(prefix + n, sel, o) if n in ("hbc", "muI") else None
        ue = film_exact(n, np.clip(zz, 0, h), h, math.radians(th), rho, par) if pp is None else film_ref_p(n, zz, pp, h, math.radians(th), rho, par)
        dz = np.diff(np.unique(np.round(zz, 9)))[0]
        film = zz <= (0.5*h if n == "muI" else h - 0.5*dz)    # mu(I): the lower half (above it p' -> 0, tau/p' -> mu_2)
        err = float(np.max(abs(uu[film] - ue[film]))/np.max(ue[film]))
        kk = np.where(film)[0][-1]
        us_m = float(uu[kk]); us_e = float(ue[kk])
        # previous output: steadiness
        t0, d0 = F[-2]
        u0 = d0["velocity"][sel, 0][o]
        rows[n] = dict(t=round(t, 2), max_err_rel=round(err, 4), u_top_of_range_model=round(us_m, 4), u_top_of_range_exact=round(us_e, 4),
                       change_last_output=round(float(np.max(abs(uu - u0)[film])/np.max(ue)), 5))
        if pp is None:
            zf = np.linspace(0, h, 300)
            ax.plot(film_exact(n, zf, h, math.radians(th), rho, par), zf*1000, color=GREY, ls="--", lw=1.2, label="exact")
        else:
            sh = zz <= h - 2*dz
            ax.plot(ue[sh], zz[sh]*1000, color=GREY, ls="--", lw=1.2, label="exact, solved p")
        ax.plot(uu[zz < h - 1e-9], zz[zz < h - 1e-9]*1000, "o", color=BLUE, ms=3, label="CFD W 130")
        if n in ("hb", "hbc"):
            G = rho*9.81*math.sin(math.radians(th)) if n == "hb" else rho*9.81*(math.sin(math.radians(th)) - 0.1*math.cos(math.radians(th)))
            ty = par[0] + (0.1*rho*9.81*math.cos(math.radians(th))*FILM_DELTA[n] if n == "hbc" else 0.0)
            ax.axhline((h - ty/G)*1000, color=ORANGE, lw=1.0, ls=":", label="yield surface")
        ax.set_title(lab, fontsize=10); ax.set_xlabel("u (m/s)"); ax.set_ylabel("z (mm)")
    axs[0].legend(loc="upper left")
    fig.tight_layout(); fig.savefig(os.path.join(out, "06_cfd_films.png"))
    return rows

def deposit(d, thr=0.0):
    """surface z(x) of phase 1 from phi along the columns of the points"""
    P = d["points"]; x = np.round(P[:, 0], 9); z = P[:, 2]; ph = d["phi"]
    xs = np.unique(x); hs = np.zeros(len(xs))
    for k, xv in enumerate(xs):
        s = x == xv; zz = z[s]; pp = ph[s]; o = np.argsort(zz); zz, pp = zz[o], pp[o]
        w = np.where(pp > thr)[0]
        if len(w) == 0 or w[0] > 0:
            continue
        i = w[-1]
        hs[k] = zz[i] + (zz[i + 1] - zz[i])*pp[i]/(pp[i] - pp[i + 1]) if i + 1 < len(zz) else zz[i]
    return xs, hs

def collapse(prefix, out):
    plt = style()
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.4))
    rows = {}
    L0 = 0.1
    for k, (n, a, col) in enumerate((("a05", 0.5, LBLUE), ("a20", 2.0, BLUE), ("a05_fine", 0.5, ORANGE))):
        if not os.path.isdir(prefix + n):
            continue
        F = cfd_frames(prefix + n)
        # front: last point with h > 10 grain diameters (5 mm); the thin tail below it (h > dx/2)
        # still creeps forward with refinement and is reported separately
        hf = 10.0*0.0005
        ts, xf, xt, um = [], [], [], []
        for t, d in F:
            xs, hs = deposit(d)
            dxg = xs[1] - xs[0]
            ts.append(t); xf.append(xs[np.where(hs > hf)[0][-1]]); xt.append(xs[np.where(hs > 0.5*dxg)[0][-1]])
            um.append(np.max(np.linalg.norm(d["velocity"], axis=1)[d["phi"] > 0]))
        xs, hs = deposit(F[-1][1])
        dL = (xf[-1] - L0)/L0
        lube = 1.2*a if a < 2.3 else 1.9*a**(2.0/3.0)
        rows[n] = dict(a=a, t_end=round(F[-1][0], 2), runout_dL_L0=round(float(dL), 3), tail_dL_L0=round(float((xt[-1] - L0)/L0), 3),
                       lube=round(lube, 3), lajeunesse=round(a, 3),
                       H_end_over_H0=round(float(hs[0]/(a*L0)), 3), umax_end=round(float(um[-1]), 4))
        lab = f"a = {a}" + (", dx 1.25 mm" if n.endswith("fine") else "")
        ax[0].plot(xs/L0, hs/L0, color=col, lw=1.6, label=lab)
        ax[0].plot([0, 1, 1], [a, a, 0], color=col, lw=0.8, ls=":")
        ax[1].plot(ts, (np.array(xf) - L0)/L0, color=col, lw=1.6, label=lab)
        ax[1].axhline(lube, color=col, lw=1.0, ls="--")
    ax[0].set_xlabel("x/L0"); ax[0].set_ylabel("z/L0"); ax[0].legend(); ax[0].set_aspect("equal")
    ax[1].set_xlabel("t (s)"); ax[1].set_ylabel(r"$(L-L_0)/L_0$  (dashed: Lube et al. 2005)"); ax[1].legend()
    fig.tight_layout(); fig.savefig(os.path.join(out, "07_cfd_collapse.png"))
    return rows

def cfdarrest(run, run_sflow, out):
    plt = style()
    g, rho, ty = 9.81, 1000.0, 89.0
    V = 0.5*0.3; a = math.sqrt(2*ty/(rho*g)); xf_ex = (3*V/(2*a))**(2.0/3.0)
    F = cfd_frames(run)
    ts, xf, um = [], [], []
    for t, d in F:
        xs, hs = deposit(d)
        ts.append(t); xf.append(xs[np.where(hs > 0.0025)[0][-1]])
        um.append(np.max(np.linalg.norm(d["velocity"], axis=1)[d["phi"] > 0]))
    xs, hs = deposit(F[-1][1])
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
    ax[0].plot(xs, hs, color=BLUE, lw=1.6, label="CFD W 130 1")
    res = dict(cfd=dict(t_end=round(F[-1][0], 2), front=round(float(xf[-1]), 4), h_wall=round(float(hs[0]), 4),
                        area=round(float(np.trapezoid(hs, xs)), 5), umax_end=round(float(um[-1]), 5)))
    if run_sflow and os.path.isdir(run_sflow):
        t2, d2 = frames(run_sflow)[-1]
        x2, h2, u2 = line(d2)
        ax[0].plot(x2, h2, color=GREEN, lw=1.4, label="SFLOW W 120 5 (case 02)")
        res["sflow"] = dict(t_end=round(t2, 1), front=round(float(front(x2, h2, 1e-3)), 4), h_wall=round(float(h2[0]), 4))
    xx = np.linspace(0, xf_ex, 300)
    ax[0].plot(xx, a*np.sqrt(xf_ex - xx), color=GREY, ls="--", lw=1.2, label=r"$h^2 = 2\tau_y(x_f-x)/\rho g$")
    ax[0].set_xlim(0, 2); ax[0].set_xlabel("x (m)"); ax[0].set_ylabel("h (m)"); ax[0].legend()
    ax[1].semilogy(ts, np.maximum(um, 1e-8), color=BLUE, lw=1.6); ax[1].set_xlabel("t (s)"); ax[1].set_ylabel("max |u| in phase 1 (m/s)")
    fig.tight_layout(); fig.savefig(os.path.join(out, "08_cfd_arrest.png"))
    res["exact"] = dict(front=round(xf_ex, 4), h_wall=round(a*math.sqrt(xf_ex), 4))
    return res

def huppert(runs, out):
    # Newtonian viscous gravity current, 2D, constant area A: x_N = 1.411 (g A^3 t/(3 nu))^(1/5)
    # (Huppert 1982); checks the no-slip wall of the implicit diffusion (D 22)
    plt = style()
    fig, ax = plt.subplots(1, 1, figsize=(5.2, 3.4))
    nu, g = 0.05, 9.81
    rows = {}
    tt = np.linspace(0.05, 10.0, 200)
    for k, run in enumerate(runs):
        F = cfd_frames(run)
        xs, hs = deposit(F[0][1]); A = np.trapezoid(hs, xs)
        ts, xf, area = [], [], []
        for t, d in F:
            xs, hs = deposit(d)
            ts.append(t); xf.append(xs[np.where(hs > 0.5*(xs[1] - xs[0]))[0][-1]]); area.append(np.trapezoid(hs, xs))
        ts = np.array(ts); xf = np.array(xf)
        xh = 1.411*(g*A**3*ts/(3.0*nu))**0.2
        m = (ts >= 1.0) & (xf < xs[-1] - 0.01)
        name = os.path.basename(run.rstrip("/"))
        rows[name] = dict(A=round(float(A), 5), x_front_end=round(float(xf[-1]), 3), x_huppert_end=round(float(xh[-1]), 3),
                          max_rel_dev_t_ge_1=round(float(np.max(np.abs(xf[m]/xh[m] - 1.0))), 3) if m.any() else None,
                          area_change=round(float(area[-1]/area[0] - 1.0), 4))
        ax.plot(ts, xf, color=(BLUE, ORANGE, LBLUE)[k % 3], lw=1.6, label=name)
        if k == 0:
            ax.plot(tt, 1.411*(g*A**3*tt/(3.0*nu))**0.2, color="k", lw=1.0, ls="--", label="Huppert (1982)")
    ax.set_xlabel("t (s)"); ax.set_ylabel("front x (m)"); ax.legend()
    fig.tight_layout(); fig.savefig(os.path.join(out, "09_viscous_current.png"))
    return rows

def semi_implicit(prefix_si, prefix_ex, out):
    # case 10: the semi-implicit SFLOW (A 210 11) against the explicit runs of the same cases and the references;
    # <prefix>mangeney, <prefix>hb_core
    plt = style()
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
    g = 9.81; th = math.radians(35); gx, gn = g*math.sin(th), g*math.cos(th); tand = math.tan(math.radians(25))
    rows = {}
    for k, (pre, ls, lab) in enumerate(((prefix_si, "-", "semi-implicit A 210 11"), (prefix_ex, ":", "explicit RK3"))):
        if pre is None or not os.path.isdir(pre + "mangeney"):
            continue
        for t, d in frames(pre + "mangeney"):
            if abs(t - 3.0) > 0.06:
                continue
            x, h, u = line(d)
            ax[0].plot(x, h, color=(BLUE, ORANGE)[k], lw=1.6, ls=ls, label=lab)
            if k == 0:
                he, ue, xb = mangeney_exact(x, t, 1.0, gn, gx, tand)
                ax[0].plot(x, he, color=GREY, lw=1.0, ls="--", label="Mangeney et al. (2000)")
        F = frames(pre + "hb_core")
        x, h, u = line(F[-1][1])
        ax[1].plot(x, h, color=(BLUE, ORANGE)[k], lw=1.6, ls=ls, label=lab)
        rows[lab] = dict(arrest_front=round(float(front(x, h, 1e-3)), 3), arrest_h_wall=round(float(h[0]), 4))
    V = 0.5*0.3; a = math.sqrt(2*89.0/(1000.0*g)); xf = (3*V/(2*a))**(2.0/3.0)
    xx = np.linspace(0, xf, 300)
    ax[1].plot(xx, a*np.sqrt(xf - xx), color=GREY, lw=1.0, ls="--", label="fully yielded")
    ax[0].set_xlim(-12, 30); ax[0].set_xlabel("x along slope (m)"); ax[0].set_ylabel("h (m)"); ax[0].set_title("01 Coulomb dam break, t = 3 s", fontsize=10); ax[0].legend()
    ax[1].set_xlim(0, 2); ax[1].set_xlabel("x (m)"); ax[1].set_ylabel("h (m)"); ax[1].set_title("02 Carbopol deposit, t = 60 s", fontsize=10); ax[1].legend()
    fig.tight_layout(); fig.savefig(os.path.join(out, "10_semi_implicit.png"))
    return rows

def lagree(prefix, out):
    # 07: the 2D benchmark of Lagree, Staron & Popinet (2011): the same mu(I) parameters (mu_s 0.32, dmu 0.28, I0 0.4,
    # d = L0/82) and a slip wall behind the column; their 2D Navier-Stokes mu(I) runs give (L - L0)/L0 = 2.2 a and
    # H/L0 = a (a < 0.5), 0.67 a^0.4 (0.5 < a < 6); contact dynamics (Staron & Hinch 2005) 2.5 a
    plt = style()
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.4))
    rows = {}
    L0 = 0.1
    for k, (n, a, col) in enumerate((("lagree_a05", 0.5, LBLUE), ("lagree_a142", 1.42, BLUE))):
        if not os.path.isdir(prefix + n):
            continue
        F = cfd_frames(prefix + n)
        ts, xf, xt = [], [], []
        for t, d in F:
            xs, hs = deposit(d)
            dxg = xs[1] - xs[0]
            ts.append(t); xt.append(xs[np.where(hs > 0.5*dxg)[0][-1]]); xf.append(xs[np.where(hs > 0.005)[0][-1]])
        xs, hs = deposit(F[-1][1])
        hl = a if a <= 0.5 else 0.67*a**0.4
        rows[n] = dict(a=a, t_end=round(F[-1][0], 2), runout_tail=round(float((xt[-1] - L0)/L0), 3),
                       runout_h5mm=round(float((xf[-1] - L0)/L0), 3), max_runout_tail=round(float((max(xt) - L0)/L0), 3),
                       NS_2p2a=round(2.2*a, 3), CD_2p5a=round(2.5*a, 3), H_L0=round(float(hs.max()/L0), 3), H_L0_lagree=round(hl, 3))
        ax[0].plot(xs/L0, hs/L0, color=col, lw=1.6, label=f"a = {a}")
        ax[0].plot([0, 1, 1], [a, a, 0], color=col, lw=0.8, ls=":")
        ax[1].plot(ts, (np.array(xt) - L0)/L0, color=col, lw=1.6, label=f"a = {a}")
        ax[1].axhline(2.2*a, color=col, lw=1.0, ls="--")
        ax[1].axhline(2.5*a, color=col, lw=1.0, ls=":")
    ax[0].set_xlabel("x/L0"); ax[0].set_ylabel("z/L0"); ax[0].legend(); ax[0].set_aspect("equal")
    ax[1].set_xlabel("t (s)"); ax[1].set_ylabel(r"$(L-L_0)/L_0$ (--: 2.2a NS, $\cdot\cdot$: 2.5a CD)"); ax[1].legend()
    fig.tight_layout(); fig.savefig(os.path.join(out, "07_cfd_collapse_lagree.png"))
    return rows

# ---------------------------------------------------------------------------------------------
# 20 - 22: the test cases of Fornes et al. (IAHR 2017, MekIT 2017, Geofaredagen 2017)

def _heights3d(d):
    """free-surface height of every vertical node column (the lowest zero crossing of phi), 3D"""
    X = d["points"]; phi = d["phi"]
    key = np.round(X[:, :2], 6); zs = X[:, 2]
    cols = {}
    for n, (a, b) in enumerate(map(tuple, key)):
        cols.setdefault((a, b), []).append(n)
    out = []
    for (a, b), idx in cols.items():
        idx = np.array(idx); o = np.argsort(zs[idx]); z = zs[idx][o]; f = phi[idx][o]
        h = 0.0
        if f[0] > 0:
            for k in range(len(z) - 1):
                if f[k] > 0 and f[k + 1] <= 0:
                    h = z[k] + (z[k + 1] - z[k])*f[k]/(f[k] - f[k + 1])
                    break     # the first crossing from the bed: not the level set extended into a lifted body
        out.append((a, b, h))
    return np.array(out)

def quickness(runs, run_sflow, out):
    # 20: quickness test of Thakur & Degago (2012), Heimdal clay c_ur = 0.2 kPa, HB parameters of Fornes et al. (2017):
    # a cylinder H 120 mm, D 100 mm released on a smooth plane; measured deposit H 40 mm, D 237 mm. Quarter domain
    # with two symmetry planes. Reference also the static fully yielded shallow pile h^2 = 2 tau_y (R - r)/(rho g)
    plt = style()
    g, rho, ty = 9.81, 2000.0, 200.0
    V = math.pi*0.05**2*0.12
    a = math.sqrt(2*ty/(rho*g)); R = (V/(2*math.pi*a*4.0/15.0))**0.4
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.4))
    rows = {}
    cols = (BLUE, LBLUE, GREEN, "#8E44AD", GREY)
    for k, (run, lab) in enumerate(runs):
        col = cols[k % len(cols)]
        if run is None or not os.path.isdir(run):
            continue
        F = cfd_frames(run)
        ts, Hs, Ds = [], [], []
        for t, d in F:
            H = _heights3d(d); w = H[:, 2] > 0.001; sx = abs(H[:, 1]) < 1e-6
            ts.append(t); Hs.append(H[:, 2].max()); Ds.append(2*H[sx & w, 0].max())
        H = _heights3d(F[-1][1]); sx = abs(H[:, 1]) < 1e-6; o = np.argsort(H[sx, 0])
        dg = abs(H[:, 0] - H[:, 1]) < 1e-6; w = H[:, 2] > 0.001
        rows[lab] = dict(t_end=round(F[-1][0], 2), H_mm=round(1e3*Hs[-1], 1), D_x_mm=round(1e3*Ds[-1], 0),
                         D_diag_mm=round(2e3*np.hypot(H[dg & w, 0], H[dg & w, 1]).max(), 0))
        ax[0].plot(1e3*H[sx, 0][o], 1e3*H[sx, 2][o], color=col, lw=1.6, label=lab)
        ax[1].plot(ts, 1e3*np.array(Ds), color=col, lw=1.6, label="D " + lab)
        ax[1].plot(ts, 1e3*np.array(Hs), color=col, lw=1.2, ls="--")
    if run_sflow is not None and os.path.isdir(run_sflow):
        F = frames(run_sflow); ts, Hs, Ds = [], [], []
        for t, d in F:
            P = d["points"]; h = d["waterlevel"]; sx = abs(P[:, 1] - P[:, 1].min()) < 1e-9
            ts.append(t); Hs.append(h.max()); Ds.append(2*P[sx & (h > 1e-3), 0].max())
        P = F[-1][1]["points"]; h = F[-1][1]["waterlevel"]; sx = abs(P[:, 1] - P[:, 1].min()) < 1e-9; o = np.argsort(P[sx, 0])
        rows["SFLOW W 120 5"] = dict(t_end=round(F[-1][0], 2), H_mm=round(1e3*Hs[-1], 1), D_x_mm=round(1e3*Ds[-1], 0))
        ax[0].plot(1e3*P[sx, 0][o], 1e3*h[sx][o], color=ORANGE, lw=1.4, label="SFLOW W 120 5 (depth-averaged)")
        ax[1].plot(ts, 1e3*np.array(Ds), color=ORANGE, lw=1.4, label="D SFLOW")
    rr = np.linspace(0, R, 200)
    ax[0].plot(1e3*rr, 1e3*a*np.sqrt(R - rr), color=GREY, ls="--", lw=1.0, label="static fully yielded")
    ax[0].plot([0, 118.5], [40, 40], color=GREEN, lw=1.0, ls=":"); ax[0].plot([118.5], [0], "v", color=GREEN, ms=7, label="measured H 40 mm, D/2 118.5 mm")
    ax[0].set_xlim(0, 160); ax[0].set_xlabel("r (mm)"); ax[0].set_ylabel("h (mm)"); ax[0].legend(fontsize=8)
    ax[1].axhline(237, color=GREEN, lw=1.0, ls=":"); ax[1].axhline(40, color=GREEN, lw=1.0, ls=":")
    ax[1].set_xlabel("t (s)"); ax[1].set_ylabel("D, H (mm) (--: H)"); ax[1].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(out, "20_quickness_test.png"))
    rows["measured"] = dict(H_mm=40, D_mm=237); rows["static fully yielded"] = dict(H_mm=round(1e3*a*math.sqrt(R), 1), D_mm=round(2e3*R, 0))
    return rows

HUANG_TOE = ((0.1, 0.370), (0.2, 0.413), (0.3, 0.475), (0.4, 0.507), (0.5, 0.51))   # digitised, h > 2 mm

def huang(runs, data, out):
    # 21: dry sand dam break of Huang et al. (2011) in the set-up of Fornes (Geofaredagen 2017): sand 0.32 x 0.20 m in
    # a box 0.95 x 0.35 m, rho 1458, rho_s 2700, d 0.24 mm; mu(I) mu_s 0.7, mu_2 1.0, I0 0.3; HBC phi 35, c 0.1 Pa,
    # K 15, n 0.35
    plt = style()
    D = np.loadtxt(data)
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.4))
    ax[0].plot(D[:, 0], D[:, 1], ".", color="#BBBBBB", ms=2, label="Huang et al. (2011), t = 0.1-0.5 s")
    rows = {}
    cols = (BLUE, ORANGE, GREEN)
    for k, run in enumerate(runs):
        if not os.path.isdir(run):
            continue
        name = os.path.basename(run.rstrip("/"))
        name = {"muI": "mu(I) W 130 3", "hbc": "HBC W 130 2", "gate": "mu(I), mu_max 1e3, gate lifted at 1 m/s", "mx": "mu(I), mu_max 1e3"}.get(name.split("_")[-1], name)
        F = cfd_frames(run); ts, toe = [], []
        for t, d in F:
            xs, hs = deposit(d); ts.append(t); toe.append(xs[np.where(hs > 0.002)[0][-1]])
            if any(abs(t - tt) < 0.026 for tt in (0.1, 0.2, 0.3, 0.4, 0.5)):
                ax[0].plot(xs, hs, color=cols[k], lw=1.2, label=name if abs(t - 0.1) < 0.026 else None)
        ts = np.array(ts); toe = np.array(toe)
        rows[name] = {("toe_%.1f" % tt): round(float(np.interp(tt, ts, toe)), 3) for tt, _ in HUANG_TOE}
        rows[name]["toe_end"] = round(float(toe[-1]), 3); rows[name]["t_end"] = round(float(ts[-1]), 2)
        ax[1].plot(ts, toe, color=cols[k], lw=1.6, label=name)
    ax[1].plot([t for t, _ in HUANG_TOE], [x for _, x in HUANG_TOE], "o", color=GREY, label="Huang et al. (2011)")
    ax[0].set_xlim(0, 0.95); ax[0].set_ylim(0, 0.25); ax[0].set_xlabel("x (m)"); ax[0].set_ylabel("z (m)"); ax[0].legend(fontsize=8)
    ax[1].set_xlabel("t (s)"); ax[1].set_ylabel("toe h > 2 mm (m)"); ax[1].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(out, "21_huang_sand.png"))
    rows["measured"] = {("toe_%.1f" % t): x for t, x in HUANG_TOE}
    return rows

def moriguchi(runs, out, width=0.3, exp=(564.0, 710.0)):
    # 22: granular flow down a 45 deg flume onto a wall (Moriguchi et al. 2009) in the 2D set-up of Fornes et al.
    # (MekIT 2017): sand 0.5 x 0.3 m, the wall 1.8 m downstream, HBC phi 30, c 1 Pa, K 1, n 1, rho 1379. Two walls:
    # the end of the domain (infinitely high, as Fornes et al.) and the 0.3 m obstacle of the experiment (S 10, the flow
    # can pass over it). The wall force is the pressure on the upstream face (gauge: minus the pressure at the top of
    # the domain) per unit width; experimental peaks of the trials at 45 deg 564 - 710 N/m (as quoted by Cuomo et al.)
    plt = style()
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.4))
    rows = {}
    for k, (run, xface, zmax, lab) in enumerate(runs):
        if not os.path.isdir(run):
            continue
        files = sorted(glob.glob(os.path.join(run, "REEF3D_CFD_VTU", "*.vtu")))
        byn = {}
        for f in files:
            m = re.search(r"-(\d+)-(\d+)\.vtu$", f)
            if m:
                byn.setdefault(int(m.group(1)), []).append(f)
        ts, Fw, snaps = [], [], []
        for n in sorted(byn):
            parts = [read_vtk(f) for f in byn[n]]
            head = open(byn[n][0], "rb").read(4000).decode("latin-1")
            t = float(re.search(r'Name="TimeValue"[^>]*>\s*([-0-9.eE+]+)', head).group(1))
            P = np.concatenate([q["points"] for q in parts]); p = np.concatenate([q["pressure"] for q in parts])
            phi = np.concatenate([q["phi"] for q in parts])
            if np.abs(np.concatenate([q["velocity"] for q in parts])).max() > 20.0:
                break     # the run blew up (EMERGENCY STOP N 61): a free fall along 2.3 m at 45 deg gives < 6 m/s
            xs = np.unique(np.round(P[:, 0], 6))
            xw = xs[-1] if xface is None else xs[xs < xface - 1e-6][-1]
            yw = P[:, 1].min()
            s = (abs(P[:, 0] - xw) < 1e-6) & (abs(P[:, 1] - yw) < 1e-9)
            z = P[s, 2]; o = np.argsort(z); z = z[o]; pw = p[s][o] - p[s][o][-1]
            sel = z <= (zmax if zmax is not None else 1e9)
            ts.append(t); Fw.append(np.trapezoid(np.maximum(pw[sel], 0.0), z[sel]))
            if any(abs(t - tt) < 0.011 for tt in (0.4, 0.8, 1.2, 1.6, 2.0)) and k == len(runs) - 1:
                snaps.append((t, {"points": P, "phi": phi}))
        ts = np.array(ts); Fw = np.array(Fw)
        ax[0].plot(ts, Fw, color=(BLUE, ORANGE, GREEN, GREY, LBLUE)[k % 5], lw=1.6, label=lab)
        i = np.argmax(Fw)
        rows[lab] = dict(peak_N_per_m=round(float(Fw[i]), 0), t_peak=round(float(ts[i]), 2), F_end_N_per_m=round(float(Fw[-1]), 0),
                         t_end=round(float(ts[-1]), 2),
                         F_at_1p2s=round(float(np.interp(1.2, ts, Fw)), 0),
                         t_contact=round(float(ts[np.argmax(Fw > 20.0)]), 2) if np.any(Fw > 20.0) else None,
                         peak_N_width_0p3=round(float(Fw[i])*width, 0))
    ax[0].axhspan(exp[0], exp[1], color=GREEN, alpha=0.15, label="experimental peaks, 45 deg")
    ax[0].set_xlabel("t (s)"); ax[0].set_ylabel("wall force (N/m)"); ax[0].legend(fontsize=8)
    for t, d in snaps:
        xs, hs = deposit(d)
        ax[1].plot(xs, hs, lw=1.2, label="t = %.1f s" % t)
    ax[1].set_xlabel("x along slope (m)"); ax[1].set_ylabel("depth normal to the slope (m)"); ax[1].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(out, "22_moriguchi_impact.png"))
    rows["experiment"] = dict(peak_N_per_m=list(exp), arrival_s="0.7 - 0.8 (as quoted)")
    return rows

# ---------------------------------------------------------------------------------------------
# 23: side-wall friction W 135 (verification)

def wallfriction(runs, out, theta=25.0, h=0.05, dw=20.0):
    # periodic layer h on a free-slip bed (symmetry plane), slope theta, two side walls at distance b (W 135 dw b):
    # the depth-integrated momentum balance with the lithostatic pressure of the uniform layer gives the mean velocity
    # u = g (sin theta - tan dw cos theta h/b) t, or rest when the wall friction exceeds the driving force
    plt = style()
    g = 9.81; th = math.radians(theta)
    fig, ax = plt.subplots(1, 1, figsize=(5.5, 3.4))
    rows = {}
    for k, run in enumerate(runs):
        if not os.path.isdir(run):
            continue
        name = os.path.basename(run.rstrip("/"))
        b = {"w_free": None, "w_b05": 0.05, "w_b03": 0.03}.get(name)
        F = cfd_frames(run); ts, um = [], []
        for t, d in F:
            x = d["points"][:, 0]
            m = (d["phi"] > 0.0) & (x > x.min() + 1e-9) & (x < x.max() - 1e-9)     # the periodic end nodes are written as 0
            ts.append(t); um.append(float(np.mean(d["velocity"][m, 0])))
        a = g*math.sin(th) - (0.0 if b is None else math.tan(math.radians(dw))*g*math.cos(th)*h/b)
        ex = np.maximum(a, 0.0)*np.array(ts)
        lab = "no wall friction" if b is None else "b = %.2f m" % b
        col = (BLUE, ORANGE, GREEN)[k % 3]
        ax.plot(ts, um, color=col, lw=1.6, label=lab + " (CFD)")
        ax.plot(ts, ex, color=col, lw=1.0, ls="--")
        rows[lab] = dict(t=round(ts[-1], 2), u_cfd=round(um[-1], 4), u_exact=round(float(ex[-1]), 4),
                         err=round(abs(um[-1] - ex[-1])/max(ex[-1], 1e-9), 4) if ex[-1] > 0 else None)
    ax.set_xlabel("t (s)"); ax.set_ylabel("mean velocity of the layer (m/s)"); ax.legend(fontsize=8)
    ax.set_title(r"dashed: $g\,(\sin\theta - \tan\delta_w \cos\theta\; h/b)\, t$", fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(out, "23_side_wall_friction.png"))
    return rows

# ---------------------------------------------------------------------------------------------
# 24: 3D periodic mu(I) film on one rank and split over two ranks across the periodic boundaries (x or y)

def periodic3d(runs, out, t=2.0):
    plt = style()
    n, th, h, rho, par, lab = FILMS[2]
    fig, ax = plt.subplots(1, 1, figsize=(5.0, 3.6))
    rows = {}; ref = None
    for k, run in enumerate(runs):
        if not os.path.isdir(run):
            continue
        F = cfd_frames(run)
        tt, d = min(F, key=lambda f: abs(f[0] - t))
        P = d["points"]; U = d["velocity"]
        xs = np.unique(np.round(P[:, 0], 7)); ys = np.unique(np.round(P[:, 1], 7))
        inn = (P[:, 0] > xs[0] + 1e-8) & (P[:, 0] < xs[-1] - 1e-8) & (P[:, 1] > ys[0] + 1e-8) & (P[:, 1] < ys[-1] - 1e-8)
        zs = np.unique(np.round(P[:, 2], 7)); pr = []; sp = 0.0
        for z in zs:
            s = inn & (abs(P[:, 2] - z) < 1e-8)
            pr.append(U[s, 0].mean()); sp = max(sp, float(U[s, 0].max() - U[s, 0].min()))
        pr = np.array(pr)
        if ref is None:
            ref = pr
        low = zs <= 0.5*h
        name = os.path.basename(run.rstrip("/"))
        rows[name] = dict(t=round(tt, 2), u_mid=round(float(np.interp(0.5*h, zs, pr)), 5), spread_horizontal=float("%.2e" % sp),
                          v_max=float("%.1e" % abs(U[inn, 1]).max()), w_max=float("%.1e" % abs(U[inn, 2]).max()),
                          max_diff_to_first_lower_half=float("%.1e" % np.abs(pr - ref)[low].max()),
                          max_diff_to_first=float("%.1e" % np.abs(pr - ref).max()))
        ax.plot(pr[zs < h - 1e-9], 1e3*zs[zs < h - 1e-9], ("o", "x", "+")[k % 3], color=(BLUE, ORANGE, GREEN)[k % 3], ms=5,
                label={"p_r1": "1 rank", "p_r2x": "2 ranks, split in x", "p_r2y": "2 ranks, split in y"}.get(name, name))
    ax.set_xlabel("u, mean over x and y (m/s)"); ax.set_ylabel("z (mm)"); ax.legend(fontsize=8)
    ax.set_title("3D periodic mu(I) film, t = %.0f s" % t, fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(out, "24_cfd_3d_periodic.png"))
    return rows

# ---------------------------------------------------------------------------------------------
# 25: 3D axisymmetric collapse of a granular column (Lube, Huppert, Sparks & Hallworth 2004), quarter domain

def axicollapse(runs, out, r0=0.05):
    # Lube et al. (2004): (r - r0)/r0 = 1.24 a (a < 1.7), 1.6 a^0.5 (a > 1.7); h/r0 = a (a < 1.7, undisturbed top),
    # 0.88 a^(1/6) (1.7 < a < 10); t_stop = 3 (h0/g)^0.5; sand rho_s 2600, d 0.32 mm, angle of repose 30 deg
    plt = style()
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.4))
    rows = {}
    for k, (run, a) in enumerate(runs):
        if not os.path.isdir(run):
            continue
        F = cfd_frames(run)
        ts, rr = [], []
        for t, d in F:
            H = _heights3d(d); sx = abs(H[:, 1]) < 1e-6; w = H[:, 2] > 0.001
            ts.append(t); rr.append(H[sx & w, 0].max())
        H = _heights3d(F[-1][1]); sx = abs(H[:, 1]) < 1e-6; dg = abs(H[:, 0] - H[:, 1]) < 1e-6; w = H[:, 2] > 0.001
        r_ax = H[sx & w, 0].max(); r_dg = np.hypot(H[dg & w, 0], H[dg & w, 1]).max()
        lube_r = 1.24*a if a < 1.7 else 1.6*a**0.5
        lube_h = a if a < 1.7 else 0.88*a**(1.0/6.0)
        rr = np.array(rr); ts = np.array(ts)
        i = np.where(np.abs(rr - rr[-1]) > 0.002*r0)[0]
        t_stop = float(ts[i[-1] + 1]) if len(i) and i[-1] + 1 < len(ts) else float(ts[0])
        rows["a = %g" % a] = dict(t_end=round(float(ts[-1]), 2), runout_axis=round(float((r_ax - r0)/r0), 3), runout_diag=round(float((r_dg - r0)/r0), 3),
                                  runout_lube=round(lube_r, 3), h_r0=round(float(H[:, 2].max()/r0), 3), h_r0_lube=round(lube_h, 3),
                                  t_stop=round(t_stop, 2), t_stop_lube=round(3.0*math.sqrt(a*r0/9.81), 3))
        o = np.argsort(H[sx, 0])
        col = (BLUE, ORANGE)[k % 2]
        ax[0].plot(H[sx, 0][o]/r0, H[sx, 2][o]/r0, color=col, lw=1.6, label="a = %g, x axis" % a)
        rd = np.hypot(H[dg, 0], H[dg, 1]); od = np.argsort(rd)
        ax[0].plot(rd[od]/r0, H[dg, 2][od]/r0, color=col, lw=1.0, ls="--", label="a = %g, diagonal" % a)
        ax[0].plot([1 + lube_r], [0], "v", color=col, ms=7)
        ax[0].plot([0, 0.3], [lube_h, lube_h], color=col, lw=1.0, ls=":")
        ax[1].plot(ts/math.sqrt(a*r0/9.81), (rr - r0)/r0, color=col, lw=1.6, label="a = %g" % a)
        ax[1].axhline(lube_r, color=col, lw=1.0, ls=":")
    ax[1].axvline(3.0, color=GREY, lw=1.0, ls="--")
    ax[0].set_xlabel("r/r0"); ax[0].set_ylabel("h/r0"); ax[0].legend(fontsize=8)
    ax[1].set_xlabel(r"$t/\sqrt{h_0/g}$"); ax[1].set_ylabel(r"$(r-r_0)/r_0$  (dotted: Lube et al. 2004)"); ax[1].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(out, "25_axisymmetric_collapse.png"))
    return rows

# ---------------------------------------------------------------------------------------------
# 26: real landslides on terrain (SFLOW W 120, initial mass from F 73): the deposit against the DEM difference

def _grid_nodes(d, nx, ny, dx, key):
    """node values of an SFLOW VTP frame (pieces merged) on the (ny+1, nx+1) node grid, the cell values the mean of
    their four nodes"""
    P = d["points"]
    i = np.clip(np.rint(P[:, 0]/dx).astype(int), 0, nx); j = np.clip(np.rint(P[:, 1]/dx).astype(int), 0, ny)
    v = d[key] if d[key].ndim == 1 else np.linalg.norm(d[key][:, :2], axis=1)
    g = np.zeros((ny + 1, nx + 1)); g[j, i] = v
    return 0.25*(g[:-1, :-1] + g[1:, :-1] + g[:-1, 1:] + g[1:, 1:])

def _drop_small(mask, min_cells):
    """mask without its patches (8 neighbours) of fewer than min_cells cells (DEM noise, small fills)"""
    try:
        from scipy import ndimage
        lab, n = ndimage.label(mask, structure=np.ones((3, 3)))
        sz = ndimage.sum(mask, lab, range(1, n + 1))
        return np.isin(lab, 1 + np.where(sz >= min_cells)[0])
    except ImportError:
        return mask

def terrainslide(run, case, out, name, hthr=0.5, fthr=0.1):
    """26: real landslide on pre/post-event DEMs (slide_terrain.py): deposit volume and area, IoU and hit rate against
    the DEM deposit (post - pre > 0.5 m outside the release area, patches of at least 1000 m^2), material left in the release
    area, runout from the release centre, front over time; optional regions.json in the case dir
    {"name": [xmin, ymin, xmax, ymax], ...} (DEM coordinates) for the volumes of branches"""
    plt = style()
    R = np.load(os.path.join(case, "reference.npz"))
    dz, rel = R["dz"], R["release"]; dx = float(R["dx"]); x0, y0 = (float(v) for v in R["origin"])
    nopost = R["nopost"] if "nopost" in R.files else np.zeros_like(rel)
    ny, nx = dz.shape; A = dx*dx
    X, Y = np.meshgrid(dx*(np.arange(nx) + 0.5), dx*(np.arange(ny) + 0.5))
    xr, yr = X[rel].mean(), Y[rel].mean()
    dist = np.hypot(X - xr, Y - yr)
    ref = _drop_small((dz > hthr) & ~rel & ~nopost, int(np.ceil(1000.0/A)))
    fr = frames(run)
    ts, front, vtot, vrel = [], [], [], []
    for t, d in fr:
        h = _grid_nodes(d, nx, ny, dx, "waterlevel")
        ts.append(t); vtot.append(h.sum()*A); vrel.append(h[rel].sum()*A)
        m = (h > fthr) & ~rel
        front.append(dist[m].max() if m.any() else 0.0)
    t, d = fr[-1]
    h = _grid_nodes(d, nx, ny, dx, "waterlevel"); U = _grid_nodes(d, nx, ny, dx, "velocity")
    mod = (h > hthr) & ~rel
    inter, union = (mod & ref).sum(), (mod | ref).sum()
    r = {"time": t, "frames": len(fr), "volume_initial_m3": vtot[0], "volume_final_m3": vtot[-1],
         "model_deposit_volume_outside_release_m3": float(h[~rel].sum()*A),
         "model_deposit_volume_h>thr_m3": float(h[mod].sum()*A),
         "dem_deposit_volume_m3": float(dz[ref].sum()*A),
         "model_deposit_area_km2": float(mod.sum()*A/1e6), "dem_deposit_area_km2": float(ref.sum()*A/1e6),
         "IoU": float(inter/max(union, 1)), "hit_rate": float(inter/max(ref.sum(), 1)),
         "false_alarm": float((mod & ~ref).sum()/max(mod.sum(), 1)),
         "model_left_in_release_m3": float(h[rel].sum()*A),
         "model_runout_m": float(dist[(h > fthr) & ~rel].max()) if ((h > fthr) & ~rel).any() else 0.0,
         "dem_runout_m": float(dist[ref].max()) if ref.any() else 0.0,
         "max_speed_final_m/s": float(U.max()), "max_depth_final_m": float(h.max())}
    regs = os.path.join(case, "regions.json")
    if os.path.exists(regs):
        for k, b in json.load(open(regs)).items():
            box = (X + x0 >= b[0]) & (X + x0 <= b[2]) & (Y + y0 >= b[1]) & (Y + y0 <= b[3]) & ~rel
            r["region_" + k] = {"model_m3": float(h[box].sum()*A), "dem_m3": float(dz[box & ref].sum()*A)}
    ext = [x0, x0 + nx*dx, y0, y0 + ny*dx]
    fig = plt.figure(figsize=(13, 8.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 0.8])
    for k, (fld, ttl) in enumerate(((np.where(ref, dz, np.nan), "DEM deposit (post - pre), m"),
                                     (np.where(h > 0.05, h, np.nan), "model depth at t = %.0f s, m" % t))):
        ax = fig.add_subplot(gs[0, k])
        ax.imshow(np.where(rel, 1.0, np.nan), origin="lower", extent=ext, cmap="Greys", vmin=0, vmax=3, alpha=0.6)
        im = ax.imshow(fld, origin="lower", extent=ext, cmap="viridis", vmin=0, vmax=max(5.0, np.nanpercentile(dz[ref], 98) if ref.any() else 5))
        ax.contour(X + x0, Y + y0, ref.astype(float), levels=[0.5], colors="#E07B39", linewidths=0.8)
        if k == 1:
            ax.contour(X + x0, Y + y0, mod.astype(float), levels=[0.5], colors="white", linewidths=0.6)
        ax.set_title(ttl, fontsize=10); ax.set_aspect("equal"); fig.colorbar(im, ax=ax, shrink=0.75)
        ax.ticklabel_format(useOffset=False, style="plain"); ax.tick_params(labelsize=7)
    ax = fig.add_subplot(gs[1, 0])
    ax.plot(ts, front, color=BLUE, label="model front (h > %.1f m)" % fthr)
    ax.axhline(r["dem_runout_m"], color=ORANGE, ls="--", label="DEM deposit runout")
    ax.set_xlabel("t (s)"); ax.set_ylabel("distance from release centre (m)"); ax.legend(fontsize=8)
    ax = fig.add_subplot(gs[1, 1])
    lab = ["deposit\noutside", "left in\nrelease"]; mv = [r["model_deposit_volume_outside_release_m3"], r["model_left_in_release_m3"]]
    dv = [r["dem_deposit_volume_m3"], 0.0]
    for k in sorted(kk for kk in r if kk.startswith("region_")):
        lab.append(k[7:]); mv.append(r[k]["model_m3"]); dv.append(r[k]["dem_m3"])
    xx = np.arange(len(lab))
    ax.bar(xx - 0.2, np.array(dv)/1e6, 0.4, color=ORANGE, label="DEM"); ax.bar(xx + 0.2, np.array(mv)/1e6, 0.4, color=BLUE, label="model")
    ax.set_xticks(xx); ax.set_xticklabels(lab, fontsize=8); ax.set_ylabel("volume (Mm$^3$)"); ax.legend(fontsize=8)
    ax.set_title("IoU %.2f, hit rate %.2f, runout %.0f / %.0f m" % (r["IoU"], r["hit_rate"], r["model_runout_m"], r["dem_runout_m"]), fontsize=9)
    fig.tight_layout(); os.makedirs(out, exist_ok=True)
    fig.savefig(os.path.join(out, "26_%s.png" % name)); plt.close(fig)
    np.savez_compressed(os.path.join(out, "26_%s_final.npz" % name), h=h, speed=U, t=np.array(ts), front=np.array(front),
                        vtot=np.array(vtot), vrel=np.array(vrel))
    return r

def terrainseries(case, out, runs, xlabel="yield stress tau_y (Pa)"):
    """26: a series of terrainslide runs (e.g. tau_y): final depth maps against the DEM deposit outline and the
    measures against the parameter; runs = [(label, value, run dir), ...]"""
    plt = style()
    res = [(lab, val, terrainslide(run, case, out, lab)) for lab, val, run in runs]
    R = np.load(os.path.join(case, "reference.npz")); dz, rel = R["dz"], R["release"]
    dx = float(R["dx"]); x0, y0 = (float(v) for v in R["origin"]); ny, nx = dz.shape
    nopost = R["nopost"] if "nopost" in R.files else np.zeros_like(rel)
    ref = _drop_small((dz > 0.5) & ~rel & ~nopost, int(np.ceil(1000.0/(dx*dx))))
    X, Y = np.meshgrid(x0 + dx*(np.arange(nx) + 0.5), y0 + dx*(np.arange(ny) + 0.5))
    n = len(res)
    fig, axs = plt.subplots(2, (n + 1)//2, figsize=(5.6*((n + 1)//2) + 1.0, 7.6), layout="constrained")
    for ax, (lab, val, r) in zip(axs.ravel(), res):
        h = np.load(os.path.join(out, "26_%s_final.npz" % lab))["h"]
        ax.imshow(np.where(rel, 1.0, np.nan), origin="lower", extent=[x0, x0 + nx*dx, y0, y0 + ny*dx], cmap="Greys", vmin=0, vmax=3, alpha=0.6)
        im = ax.imshow(np.where(h > 0.05, h, np.nan), origin="lower", extent=[x0, x0 + nx*dx, y0, y0 + ny*dx], cmap="viridis", vmin=0, vmax=10)
        ax.contour(X, Y, ref.astype(float), levels=[0.5], colors=ORANGE, linewidths=0.8)
        ax.set_title("%s: IoU %.2f, hit %.2f, runout %.0f m" % (lab, r["IoU"], r["hit_rate"], r["model_runout_m"]), fontsize=9)
        ax.set_aspect("equal"); ax.ticklabel_format(useOffset=False, style="plain"); ax.tick_params(labelsize=7)
    fig.colorbar(im, ax=axs.ravel().tolist(), shrink=0.5, label="final depth (m)")
    os.makedirs(out, exist_ok=True); fig.savefig(os.path.join(out, "26_series_maps.png")); plt.close(fig)
    v = np.array([q[1] for q in res])
    fig, axs = plt.subplots(1, 3, figsize=(13, 3.8))
    axs[0].plot(v, [q[2]["model_runout_m"] for q in res], "o-", color=BLUE, label="model")
    axs[0].axhline(res[0][2]["dem_runout_m"], color=ORANGE, ls="--", label="DEM"); axs[0].set_ylabel("runout from release centre (m)")
    axs[1].plot(v, [q[2]["IoU"] for q in res], "o-", color=BLUE, label="IoU")
    axs[1].plot(v, [q[2]["hit_rate"] for q in res], "s-", color=LBLUE, label="hit rate")
    axs[1].plot(v, [q[2]["false_alarm"] for q in res], "^-", color=GREY, label="false alarm"); axs[1].set_ylim(0, 1)
    regs = sorted(k for k in res[0][2] if k.startswith("region_"))
    cols = [BLUE, GREEN, LBLUE, GREY]
    for c, k in zip(cols, regs):
        axs[2].plot(v, [q[2][k]["model_m3"]/1e6 for q in res], "o-", color=c, label=k[7:])
        axs[2].axhline(res[0][2][k]["dem_m3"]/1e6, color=c, ls="--", lw=0.8)
    axs[2].plot(v, [q[2]["model_left_in_release_m3"]/1e6 for q in res], "o-", color=ORANGE, label="left in release")
    axs[2].set_ylabel("volume (Mm$^3$), DEM dashed")
    axs[0].set_title("runout", fontsize=10); axs[1].set_title("footprint against the DEM deposit", fontsize=10)
    axs[2].set_title("volumes of the branches", fontsize=10)
    for ax in axs:
        ax.set_xscale("log"); ax.minorticks_off(); ax.set_xticks(v); ax.set_xticklabels(["%g" % q for q in v])
        ax.set_xlabel(xlabel); ax.legend(fontsize=8, loc="best")
    fig.tight_layout(); fig.savefig(os.path.join(out, "26_series_measures.png")); plt.close(fig)
    return {lab: r for lab, val, r in res}

def _sflow_max_cells(run, nx, ny, dx):
    """REEF3D_SFLOW_MAX (P 112) on the cell grid: max_eta, max_depth, max_speed, time_max_eta and inundation of each
    cell (the two triangles of a cell, merged pieces); NaN where never wet"""
    files = sorted(glob.glob(os.path.join(run, "REEF3D_SFLOW_MAX", "REEF3D-SFLOW-MAX-*.vtp")))
    out = {k: np.full((ny, nx), np.nan) for k in ("max_eta", "max_depth", "max_speed", "time_max_eta", "inundation")}
    for fn in files:
        d = read_vtk(fn)
        tri = d["connectivity"].reshape(-1, 3)
        c = d["points"][tri].mean(axis=1)
        i = np.clip((c[:, 0]/dx).astype(int), 0, nx - 1); j = np.clip((c[:, 1]/dx).astype(int), 0, ny - 1)
        for k in out:
            v = d[k].astype(float)
            ok = np.isfinite(v)
            a = out[k]
            cur = a[j[ok], i[ok]]
            a[j[ok], i[ok]] = np.where(np.isnan(cur), v[ok], np.maximum(cur, v[ok]))
    return out

def tsunamirunup(run, case, obs_csv, out, name, radius=100.0, slidepad=300.0):
    """27: landslide tsunami (two-layer W 150, P 112 maxima) against a runup survey: for every observed runup (csv with
    x, y in the coordinates of the DEM, flow_elev, type, location) the highest water level max_eta of the cells wet at
    any time within radius of the point (the water reaching up the slope near the mark); points inside the slide area
    (+ slidepad) are left out. Aida K = exp(mean ln(obs/model)) and kappa = exp(std ln(obs/model)) (K ~ 1 and
    kappa < 1.45: good agreement for tsunami runup) over the full runups the wave reached (model above 0.5 m)"""
    import csv
    plt = style()
    R = np.load(os.path.join(case, "reference.npz"))
    rel = R["release"]; dx = float(R["dx"]); x0, y0 = (float(v) for v in R["origin"]); ny, nx = rel.shape
    pre = R["pre"]
    M = _sflow_max_cells(run, nx, ny, dx)
    X, Y = np.meshgrid(dx*(np.arange(nx) + 0.5), dx*(np.arange(ny) + 0.5))
    xr, yr = X[rel].mean(), Y[rel].mean()
    dist_rel = np.full(rel.shape, np.inf)
    ry, rx = np.where(rel)
    rows = [r for r in csv.DictReader(open(obs_csv))]
    obs = []
    for r in rows:
        if not r["type"].startswith("Runup") or r["type"] == "Runup_questionable" or not r["flow_elev"].strip():
            continue
        if float(r["flow_elev"]) < 0.5:        # marks at or below the sea level (e.g. trimline points of a DEM)
            continue
        x, y = float(r["x"]) - x0, float(r["y"]) - y0
        if not (0 <= x < nx*dx and 0 <= y < ny*dx):
            continue
        if np.min(np.hypot(X[rel] - x, Y[rel] - y)) < slidepad:
            continue
        obs.append((x, y, float(r["flow_elev"]), r["type"], r["location"]))
    wet = np.isfinite(M["max_eta"])
    res = []
    for x, y, z, ty, loc in obs:
        i0, i1 = int((x - radius)/dx), int((x + radius)/dx) + 1
        j0, j1 = int((y - radius)/dx), int((y + radius)/dx) + 1
        sl = (slice(max(j0, 0), min(j1 + 1, ny)), slice(max(i0, 0), min(i1 + 1, nx)))
        near = (np.hypot(X[sl] - x, Y[sl] - y) <= radius) & wet[sl]
        zm = float(np.max(M["max_eta"][sl][near])) if near.any() else 0.0
        res.append((x, y, z, zm, ty, loc, float(np.hypot(x - xr, y - yr))))
    A = np.array([(q[2], q[3], q[6]) for q in res])
    full = np.array([q[4] == "Runup" for q in res])
    o, m = A[:, 0], np.maximum(A[:, 1], 0.1)
    reached = A[:, 1] > 0.5               # the wave came within the radius above 0.5 m (else not reached in the run)
    lr = np.log(o[full & reached]/m[full & reached])
    r_out = {"points": len(res), "points_runup": int(full.sum()), "points_runup_min": int((~full).sum()),
             "obs_max_m": float(o.max()), "model_max_near_obs_m": float(A[:, 1].max()),
             "model_max_eta_m": float(np.nanmax(M["max_eta"][~rel])),
             "points_not_reached": int((~reached).sum()),
             "Aida_K": float(np.exp(lr.mean())), "Aida_kappa": float(np.exp(lr.std())),
             "within_factor_1.5": float(np.mean(np.abs(lr) < np.log(1.5))),
             "runup_min_exceeded": float(np.mean(A[~full, 1] >= A[~full, 0])) if (~full).any() else None,
             "flooded_area_km2": float(np.sum((M["inundation"] == 2) & (pre > 0))*dx*dx/1e6),
             "radius_m": radius}
    locs = {}
    for q, f_ in zip(res, full):
        if f_:
            locs.setdefault(q[5], []).append((q[2], q[3]))
    r_out["by_location_max_obs_model"] = {k: [round(max(a for a, b in v), 1), round(max(b for a, b in v), 1)] for k, v in locs.items()}
    fig = plt.figure(figsize=(14, 9))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.15, 1])
    ax = fig.add_subplot(gs[:, 0])
    from matplotlib.colors import LightSource
    hs = LightSource(azdeg=315, altdeg=45).hillshade(np.nan_to_num(pre[::-1]), vert_exag=1, dx=dx, dy=dx)[::-1]
    ext = [x0, x0 + nx*dx, y0, y0 + ny*dx]
    ax.imshow(hs, origin="lower", extent=ext, cmap="gray", vmin=-0.2, vmax=1.3)
    me = np.where((M["inundation"] == 2) & (pre > 0), M["max_eta"], np.nan)
    im = ax.imshow(me, origin="lower", extent=ext, cmap="viridis", vmin=0, vmax=200)
    ax.contour(X + x0, Y + y0, rel.astype(float), levels=[0.5], colors="white", linewidths=0.8)
    sc = ax.scatter(A[:, 0]*0 + np.array([q[0] for q in res]) + x0, np.array([q[1] for q in res]) + y0, c=o, s=10,
                    cmap="viridis", vmin=0, vmax=200, edgecolors="k", linewidths=0.3)
    ax.set_title("model: highest water level on flooded land (m); dots: survey", fontsize=10)
    ax.set_aspect("equal"); ax.ticklabel_format(useOffset=False, style="plain"); ax.tick_params(labelsize=7)
    fig.colorbar(im, ax=ax, shrink=0.6, label="elevation (m)")
    ax = fig.add_subplot(gs[0, 1])
    ax.plot(A[full, 2]/1000, o[full], "o", ms=3, color=ORANGE, label="survey runup")
    ax.plot(A[~full, 2]/1000, o[~full], "^", ms=3, mfc="none", color=ORANGE, label="survey runup (minimum)")
    ax.plot(A[:, 2]/1000, A[:, 1], "s", ms=2.5, color=BLUE, label="model (wet within %.0f m)" % radius)
    ax.set_xlabel("distance from the slide (km)"); ax.set_ylabel("runup (m)"); ax.set_yscale("log"); ax.set_ylim(1, 300)
    ax.legend(fontsize=8)
    ax = fig.add_subplot(gs[1, 1])
    ax.loglog(o[full], m[full], "o", ms=3, color=BLUE, label="runup")
    ax.loglog(o[~full], m[~full], "^", ms=3, mfc="none", color=GREY, label="runup (minimum)")
    ax.plot([1, 300], [1, 300], "k-", lw=0.8); ax.plot([1, 300], [1.5, 450], "k:", lw=0.6); ax.plot([1, 300], [1/1.5, 200], "k:", lw=0.6)
    ax.set_xlim(1, 300); ax.set_ylim(1, 300); ax.set_xlabel("survey (m)"); ax.set_ylabel("model (m)")
    ax.set_title("Aida K %.2f, kappa %.2f, within x1.5: %.0f %% (%d runup points reached)" % (r_out["Aida_K"], r_out["Aida_kappa"],
                 100*r_out["within_factor_1.5"], int((full & reached).sum())), fontsize=9)
    ax.legend(fontsize=8)
    fig.tight_layout(); os.makedirs(out, exist_ok=True)
    stem = name if re.match(r"^\d+_", name) else "27_" + name       # a name with a case number keeps it (28_...)
    fig.savefig(os.path.join(out, stem + ".png")); plt.close(fig)
    np.savez_compressed(os.path.join(out, stem + "_points.npz"), pts=A, full=full)
    return r_out

def tsunamiseries(case, obs_csv, out, runs, radius=100.0):
    """27: a series of tsunamirunup runs (e.g. the basal friction): runup against the distance from the slide for every
    run with the survey, and the highest runup per location; runs = [(label, run dir), ...]"""
    plt = style()
    stem = lambda lab: lab if re.match(r"^\d+_", lab) else "27_" + lab
    res = [(lab, tsunamirunup(run, case, obs_csv, out, lab, radius),
            np.load(os.path.join(out, stem(lab) + "_points.npz"))) for lab, run in runs]
    fig, axs = plt.subplots(1, 2, figsize=(14, 4.8), gridspec_kw={"width_ratios": [1.3, 1]})
    P0 = res[0][2]; A0, f0 = P0["pts"], P0["full"]
    axs[0].plot(A0[f0, 2]/1000, A0[f0, 0], "o", ms=3.5, color=ORANGE, alpha=0.6, label="survey runup", zorder=1)
    axs[0].plot(A0[~f0, 2]/1000, A0[~f0, 0], "^", ms=3.5, mfc="none", color=ORANGE, label="survey runup (minimum)", zorder=1)
    cols = [BLUE, GREEN, LBLUE, GREY]
    for c, (lab, r, P) in zip(cols, res):
        A = P["pts"]; o = np.argsort(A[:, 2])
        axs[0].plot(A[o, 2]/1000, A[o, 1], "s", ms=2, color=c, alpha=0.8, zorder=3,
                    label="%s: K %.2f, kappa %.2f" % (lab, r["Aida_K"], r["Aida_kappa"]))
    axs[0].set_yscale("log"); axs[0].set_ylim(1, 300); axs[0].set_xlabel("distance from the slide (km)"); axs[0].set_ylabel("runup (m)")
    axs[0].legend(fontsize=8)
    locs = sorted(res[0][1]["by_location_max_obs_model"], key=lambda k: -res[0][1]["by_location_max_obs_model"][k][0])
    xx = np.arange(len(locs)); w = 0.8/(len(res) + 1)
    axs[1].bar(xx - 0.4 + w/2, [res[0][1]["by_location_max_obs_model"][k][0] for k in locs], w, color=ORANGE, label="survey")
    for n, (c, (lab, r, P)) in enumerate(zip(cols, res)):
        axs[1].bar(xx - 0.4 + w*(n + 1.5), [r["by_location_max_obs_model"][k][1] for k in locs], w, color=c, label=lab)
    axs[1].set_xticks(xx); axs[1].set_xticklabels(locs, rotation=35, ha="right", fontsize=8); axs[1].set_ylabel("highest runup (m)")
    axs[1].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(out, (re.match(r"^(\d+)_", runs[0][0]).group(1) if re.match(r"^\d+_", runs[0][0]) else "27") + "_series.png")); plt.close(fig)
    return {lab: r for lab, r, P in res}

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
    elif cmd == "films":
        r = films(sys.argv[2], sys.argv[3])
    elif cmd == "collapse":
        r = collapse(sys.argv[2], sys.argv[3])
    elif cmd == "cfdarrest":
        r = cfdarrest(sys.argv[2], sys.argv[3] if len(sys.argv) > 4 else None, sys.argv[-1])
    elif cmd == "quickness":
        # quickness <out> <sflow run|-> <cfd run>:<label> ...
        r = quickness([tuple(a.split(":", 1)) for a in sys.argv[4:]], sys.argv[3] if sys.argv[3] != "-" else None, sys.argv[2])
    elif cmd == "huang":
        r = huang(sys.argv[2:-2], sys.argv[-2], sys.argv[-1])
    elif cmd == "terrainseries":
        r = terrainseries(sys.argv[2], sys.argv[3], [(a.split(":")[0], float(a.split(":")[1]), a.split(":")[2]) for a in sys.argv[4:]])
    elif cmd == "tsunamiseries":
        r = tsunamiseries(sys.argv[2], sys.argv[3], sys.argv[4], [tuple(a.split(":", 1)) for a in sys.argv[5:]])
    elif cmd == "tsunamirunup":
        r = tsunamirunup(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], *(float(a) for a in sys.argv[7:8]))
    elif cmd == "terrainslide":
        # terrainslide <run> <case dir with reference.npz> <out> <name>
        r = terrainslide(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    elif cmd == "periodic3d":
        r = periodic3d(sys.argv[2:-1], sys.argv[-1])
    elif cmd == "axicollapse":
        # axicollapse <out> <run>:<a> ...
        r = axicollapse([(a.split(":")[0], float(a.split(":")[1])) for a in sys.argv[3:]], sys.argv[2])
    elif cmd == "wallfriction":
        r = wallfriction(sys.argv[2:-1], sys.argv[-1])
    elif cmd == "moriguchi":
        # moriguchi <out> <run>:<x of the wall face|->:<height of the wall|->:<label> ...   (-: the end of the domain)
        rr = []
        for a in sys.argv[3:]:
            run, xf, zm, lab = a.split(":", 3)
            rr.append((run, None if xf == "-" else float(xf), None if zm == "-" else float(zm), lab))
        r = moriguchi(rr, sys.argv[2])
    elif cmd == "lagree":
        r = lagree(sys.argv[2], sys.argv[3])
    elif cmd == "si":
        r = semi_implicit(sys.argv[2], sys.argv[3] if len(sys.argv) > 4 else None, sys.argv[-1])
    elif cmd == "huppert":
        r = huppert(sys.argv[2:-1], sys.argv[-1])
    print(json.dumps(r, indent=1))
