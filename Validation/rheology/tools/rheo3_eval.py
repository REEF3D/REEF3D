#!/usr/bin/env python3
# Architect: Hans Bihs
"""Evaluation of the CFD three-phase cases of Validation/rheology (D 24 dynamic diffusion, W 140 three-phase slide).

  rheo3_eval.py filmair <run D24 0> <run D24 1> <run D24 2> <run no air> <out>   11: HB film under air, kinematic vs dynamic
  rheo3_eval.py current <run> <out>                       12a: viscous gravity current of 09 as a three-phase slide
  rheo3_eval.py collapse <run base> <run D24 1> <run 3phase> <out>   12b: mu(I) collapse a = 0.5 of 07: two-phase, dynamic, three-phase
  rheo3_eval.py pelayer <run> [<run fine>] <out>          13: excess pore pressure of a layer at rest vs the series solution
  rheo3_eval.py rondon <run loose> <run dense> <out>      14: granular collapse in a viscous fluid, loose and dense (Rondon et al. 2011)
  rheo3_eval.py robbe <run 10 mm> <run 5 mm> <out>        15: granular column collapse into water (Robbe-Saule et al. 2021)

Reads the CFD VTU files (P 18 2) with ../../turbulence/tools/r3read.py and the slide volume of reef3d.log.
"""
import sys, os, re, glob, math, json
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "turbulence", "tools"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r3read import read_vtk
import rheo_eval as R

BLUE, LBLUE, ORANGE, GREEN, GREY = R.BLUE, R.LBLUE, R.ORANGE, R.GREEN, R.GREY
G = 9.81


def frames(run, keys=("points", "velocity", "phi", "pressure", "alpha_s", "p_eff", "p_hyd", "p_e", "phi_s"),
           last_only=False, every=1):
    files = sorted(glob.glob(os.path.join(run, "REEF3D_CFD_VTU", "*.vtu")))
    byn = {}
    for f in files:
        m = re.search(r"-(\d+)-(\d+)\.vtu$", f)
        if m:
            byn.setdefault(int(m.group(1)), []).append(f)
    ks = sorted(byn)
    ks = ks[-1:] if last_only else ks[::every]
    out = []
    for n in ks:
        parts = [read_vtk(f) for f in byn[n]]
        head = open(byn[n][0], "rb").read(4000).decode("latin-1")
        mt = re.search(r'Name="TimeValue"[^>]*>\s*([-0-9.eE+]+)', head)
        t = float(mt.group(1)) if mt else float(n)
        d = {k: np.concatenate([p[k] for p in parts]) for k in keys if k in parts[0]}
        out.append((t, d))
    return out


def columns(d, field):
    """per x column of the middle y plane (2D: the y = min plane; the planes on side walls average with the
    ghost cells): z sorted and the field"""
    P = d["points"]; x = np.round(P[:, 0], 9); y = np.round(P[:, 1], 9); z = P[:, 2]
    ys = np.unique(y)
    s0 = y == ys[len(ys)//2] if len(ys) > 2 else y == ys[0]
    xs = np.unique(x[s0]); out = []
    for xv in xs:
        s = s0 & (x == xv)
        zz = z[s]; ff = d[field][s]
        zz, idx = np.unique(np.round(zz, 9), return_index=True)
        out.append((xv, zz, ff[idx]))
    return out


def slide_height(d):
    """h(x) = integral of alpha_s dz over the column. The nodes on walls average the cells with the solid ghost
    cells (alpha 0 there): the nodes on the bed are replaced by the first node above, the columns on the side
    walls by the next column."""
    xs, hs = [], []
    for xv, zz, aa in columns(d, "alpha_s"):
        aa = np.clip(aa, 0, 1).copy()
        aa[0] = aa[1]
        xs.append(xv); hs.append(float(np.sum(0.5*(aa[1:] + aa[:-1])*np.diff(zz))))
    xs = np.array(xs); hs = np.array(hs)
    hs[0] = hs[1]; hs[-1] = hs[-2]
    return xs, hs


def water_surface(d):
    """water surface of each column: the highest phi = 0 crossing with alpha_s < 1/2 below it"""
    xs, es = [], []
    al = {xv: aa for xv, zz, aa in columns(d, "alpha_s")}
    for xv, zz, ph in columns(d, "phi"):
        aa = al[xv]
        w = np.where((ph > 0) & (aa < 0.5))[0]
        if len(w) == 0:
            xs.append(xv); es.append(np.nan); continue
        i = w[-1]
        e = zz[i] + (zz[i + 1] - zz[i])*ph[i]/(ph[i] - ph[i + 1]) if i + 1 < len(zz) and ph[i + 1] < 0 else zz[i]
        xs.append(xv); es.append(e)
    return np.array(xs), np.array(es)


def front(x, h, thr):
    w = np.where(h > thr)[0]
    return x[w[-1]] if len(w) else np.nan


def log_volume(run):
    v = []
    for L in open(os.path.join(run, "reef3d.log"), errors="ignore"):
        m = re.match(r"Slide volume: ([-0-9.eE+]+) \(rel. change ([-0-9.eE+]+)\), alpha_s ([-0-9.eE+]+) - ([-0-9.eE+]+), clipped ([-0-9.eE+]+)", L)
        if m:
            v.append([float(g) for g in m.groups()])
    return np.array(v)


def wall(run):
    for L in open(os.path.join(run, "reef3d.log"), errors="ignore"):
        if L.startswith("wall "):
            return int(L.split()[1])
    return None


# ---------------------------------------------------------------------------------------------------------
def filmair(r0, r1, r2, rna, out):
    plt = R.style()
    h, th = 0.1, math.radians(12)
    par = [89.0, 47.68, 0.415]
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    rows = {}
    zf = np.linspace(0, h, 300)
    ax[0].plot(R.film_exact("hb", zf, h, th, 1000.0, par), zf*1000, color=GREY, ls="--", lw=1.2, label="exact (single layer)")
    for run, lab, col, mk in ((r0, "D 24 0 (kinematic)", ORANGE, "s"), (r1, "D 24 1 (dynamic)", BLUE, "o"),
                              (r2, "D 24 2 (dynamic, harmonic)", GREEN, "^"), (rna, "no air (06), dynamic", LBLUE, "x")):
        if run is None or not os.path.isdir(run):
            continue
        F = frames(run, keys=("points", "velocity", "phi"))
        t, d = F[-1]
        P = d["points"]; x = P[:, 0]; z = P[:, 2]
        xs = np.unique(np.round(x, 9)); sel = abs(x - xs[len(xs)//2]) < 1e-9
        zz, idx = np.unique(np.round(z[sel], 9), return_index=True); uu = d["velocity"][sel, 0][idx]
        ue = R.film_exact("hb", np.clip(zz, 0, h), h, th, 1000.0, par)
        film = (zz > 1e-9) & (zz <= h - 1e-9)        # without the node on the bed (interpolated with the ghost cell)
        err = float(np.max(abs(uu[film] - ue[film]))/np.max(ue))
        t0, d0 = F[-2]
        u0 = d0["velocity"][sel, 0][idx]
        rows[lab] = dict(t=round(t, 2), max_err_rel=round(err, 4), u_surface=round(float(uu[film][-1]), 4),
                         u_surface_exact=round(float(ue[film][-1]), 4), air_max=round(float(uu.max()), 3),
                         change_last_output=float("%.2e" % (np.max(abs(uu - u0)[film])/np.max(ue))), wall_s=wall(run))
        ax[0].plot(uu[zz <= h + 1e-9], zz[zz <= h + 1e-9]*1000, mk, color=col, ms=3.5, mfc="none", label=lab)
        ax[1].plot(uu, zz*1000, "-", color=col, lw=1.4, label=lab)
    ax[0].set_xlabel("u (m/s)"); ax[0].set_ylabel("z (mm)"); ax[0].set_title("film, t = 5 s", fontsize=10); ax[0].legend(loc="upper left", fontsize=8)
    ax[1].axhline(100, color=GREY, lw=0.8, ls=":"); ax[1].set_xlabel("u (m/s)"); ax[1].set_ylabel("z (mm)")
    ax[1].set_title("film and air (no-slip lid at 120 mm)", fontsize=10); ax[1].set_xscale("symlog", linthresh=0.2)
    fig.tight_layout(); fig.savefig(os.path.join(out, "11_film_air.png"))
    return rows


def current(run, run2p_note, out):
    plt = R.style()
    nu, A = 0.05, 0.005
    F = frames(run, keys=("points", "velocity", "alpha_s"))
    ts, xf = [], []
    for t, d in F:
        x, h = slide_height(d)
        ts.append(t); xf.append(front(x, h, 0.5*0.0025))
    ts = np.array(ts); xf = np.array(xf)
    hup = lambda t: 1.411*(G*A**3*t/(3*nu))**0.2
    V = log_volume(run)
    rows = dict(front_1s=round(float(np.interp(1.0, ts, xf)), 4), front_10s=round(float(xf[-1]), 4),
                huppert_1s=round(hup(1.0), 4), huppert_10s=round(hup(ts[-1]), 4),
                max_dev_t_ge_1=round(float(np.max(abs(xf[ts >= 1]/hup(ts[ts >= 1]) - 1))), 3),
                volume_rel_change_max=float("%.2e" % np.max(abs(V[:, 1]))) if len(V) else None,
                alpha_min=float(V[:, 2].min()) if len(V) else None, alpha_max=float(V[:, 3].max()) if len(V) else None,
                wall_s=wall(run))
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.4))
    tt = np.linspace(0.3, ts[-1], 200)
    ax[0].plot(ts, xf, "o-", color=BLUE, ms=3, label="three-phase W 140 (alpha_s)")
    ax[0].plot([1.0, 10.0], [0.147, 0.225], "s", color=ORANGE, label="two-phase level set (case 09)")
    ax[0].plot(tt, hup(tt), "--", color=GREY, label="Huppert (1982)")
    ax[0].set_xlabel("t (s)"); ax[0].set_ylabel("front (m)"); ax[0].legend()
    x, h = slide_height(F[-1][1])
    ax[1].plot(x, h*1000, color=BLUE, lw=1.6, label="t = %.0f s" % F[-1][0])
    x1, h1 = slide_height(F[min(1, len(F) - 1)][1])
    ax[1].plot(x1, h1*1000, color=LBLUE, lw=1.2, label="t = %.0f s" % F[min(1, len(F) - 1)][0])
    ax[1].set_xlabel("x (m)"); ax[1].set_ylabel(r"$\int \alpha_s\,dz$ (mm)"); ax[1].legend()
    fig.tight_layout(); fig.savefig(os.path.join(out, "12a_viscous_current_3phase.png"))
    return rows


def deposit_any(d):
    if "alpha_s" in d:
        return slide_height(d)
    return R.deposit(d)


def collapse(rbase, rdyn, r3, out):
    plt = R.style()
    L0, a = 0.1, 0.5
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.4))
    rows = {}
    for run, lab, col in ((rbase, "two-phase, kinematic (07)", GREY), (rdyn, "two-phase, D 24 2", BLUE), (r3, "three-phase W 140", ORANGE)):
        if run is None or not os.path.isdir(run):
            continue
        F = frames(run, keys=("points", "velocity", "phi", "alpha_s"))
        ts, xf = [], []
        for t, d in F:
            xs, hs = deposit_any(d)
            ts.append(t); xf.append(front(xs, hs, 0.005))
        xs, hs = deposit_any(F[-1][1])
        rows[lab] = dict(t_end=round(F[-1][0], 2), runout_dL_L0=round((xf[-1] - L0)/L0, 3),
                         tail_2mm_dL_L0=round((front(xs, hs, 0.002) - L0)/L0, 3), H_end_over_H0=round(float(hs[0]/(a*L0)), 3),
                         wall_s=wall(run))
        V = log_volume(run) if os.path.exists(os.path.join(run, "reef3d.log")) else []
        if len(V):
            rows[lab]["volume_rel_change_max"] = float("%.2e" % np.max(abs(V[:, 1])))
        ax[0].plot(xs/L0, hs/L0, color=col, lw=1.6, label=lab)
        ax[1].plot(ts, (np.array(xf) - L0)/L0, color=col, lw=1.6, label=lab)
    ax[0].plot([0, 1, 1], [a, a, 0], color=GREY, lw=0.8, ls=":")
    ax[1].axhline(1.2*a, color=GREY, ls="--", lw=1.0)
    ax[0].set_xlabel("x/L0"); ax[0].set_ylabel("h/L0"); ax[0].legend(fontsize=8); ax[0].set_xlim(0, 2.2)
    ax[1].set_xlabel("t (s)"); ax[1].set_ylabel(r"$(L-L_0)/L_0$  (dashed: Lube et al. 2005)"); ax[1].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(out, "12b_collapse_3phase.png"))
    return rows


def pe_series(z, t, h, D, n=200):
    """p_e/p_e0(z=0) for p_e(z,0) = p0 (1 - z/h), p_e(h) = 0, dp/dz(0) = 0"""
    out = np.zeros_like(np.atleast_1d(z), dtype=float)
    for m in range(n):
        k = (2*m + 1)*math.pi/(2*h)
        # coefficient of the linear profile 1 - z/h on cos(k z): 2/h int_0^h (1 - z/h) cos(k z) dz = 2/(h^2 k^2)
        b = 2.0/(h*h*k*k)
        out += b*np.cos(k*np.atleast_1d(z))*math.exp(-k*k*D*t)
    return out


def pelayer(runs, out):
    plt = R.style()
    h, lam, drho = 0.1, 0.5, 1000.0
    D = 1.0e-9/(1.0e-3*1.2337e-3)
    T = 4*h*h/(math.pi**2*D)
    p0 = lam*drho*G*h
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
    rows = {}
    cols = (BLUE, ORANGE)
    zq = 0.01          # a node clear of the bed (the node on the bed is averaged with the ghost cell)
    tt = np.linspace(0, 10, 300)
    ax[0].plot(tt/T, [pe_series(zq, t, h, D)[0] for t in tt], "--", color=GREY, label="series, drained at z = h")
    for k, run in enumerate(runs):
        F = frames(run, keys=("points", "velocity", "alpha_s", "p_e", "p_eff", "p_hyd", "pressure"))
        P = F[0][1]["points"]; x = np.round(P[:, 0], 9)
        sel = (x == np.unique(x)[len(np.unique(x))//2]) & (P[:, 1] == P[:, 1].min())
        zz, idx = np.unique(np.round(P[sel, 2], 9), return_index=True)
        dz = zz[1] - zz[0]
        iq = int(np.argmin(abs(zz - zq)))
        ts = np.array([t for t, d in F]); pq = np.array([d["p_e"][sel][idx][iq] for t, d in F])
        ex = np.array([pe_series(zq, t, h, D)[0] for t in ts])*p0
        exh = np.array([pe_series(zq, t, h + dz/2, D)[0]*(h + dz/2)/h for t in ts])*p0
        m = ts > 0.1
        # level 1: p - p_hyd against (rho_s - rho_w) g (h - z) at the end (p_e decayed: p_eff + p_e)
        t, d = F[-1]
        pf = d["p_eff"][sel][idx]; pe = d["p_e"][sel][idx]
        inner = (zz > 1.5*dz) & (zz < h - 1.5*dz)
        lith = drho*G*(h - zz)
        rows[run] = dict(dz=round(dz, 4), T=round(T, 3),
                         pe_z01_T=round(float(np.interp(T, ts, pq)), 1), series_T=round(float(np.interp(T, ts, ex)), 1),
                         series_T_drained_half_cell_out=round(float(np.interp(T, ts, exh)), 1),
                         max_dev_from_series_rel_p0=round(float(np.max(abs(pq[m] - ex[m]))/p0), 4),
                         max_dev_from_shifted_series_rel_p0=round(float(np.max(abs(pq[m] - exh[m]))/p0), 4),
                         p_minus_phyd_max_err_rel=round(float(np.max(abs(pf[inner] + pe[inner] - lith[inner]))/(drho*G*h)), 4), wall_s=wall(run))
        ax[0].plot(ts/T, pq/p0, "o", color=cols[k], ms=3, mfc="none", label="CFD dz = %.1f mm" % (dz*1000))
        s = zz <= h + 1e-9
        ax[1].plot((pf[s] + pe[s])/1000, zz[s]*1000, "-", color=cols[k], lw=1.4, label="p - p_hyd, dz %.1f mm" % (dz*1000))
        t1, d1 = F[1]
        ax[1].plot(d1["p_e"][sel][idx][s]/1000, zz[s]*1000, ":", color=cols[k], lw=1.4, label="p_e, t = %.1f s" % t1)
    zf = np.linspace(0, h, 50)
    ax[1].plot(drho*G*(h - zf)/1000, zf*1000, "--", color=GREY, label=r"$(\rho_s-\rho_w) g (h-z)$")
    ax[0].set_xlabel(r"$t/T$,  $T = 4h^2/(\pi^2 D)$"); ax[0].set_ylabel(r"$p_e/(\lambda_0 \rho' g h)$ at z = 10 mm"); ax[0].legend(fontsize=8)
    ax[1].set_xlabel("pressure (kPa)"); ax[1].set_ylabel("z (mm)"); ax[1].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(out, "13_pore_pressure_layer.png"))
    return rows


def rondon(rl, rd, out):
    plt = R.style()
    L0 = 0.06
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.5))
    rows = {}
    for run, lab, col, H0 in ((rl, "loose, phi 0.55", ORANGE, 0.048), (rd, "dense, phi 0.60", BLUE, 0.042)):
        if run is None or not os.path.isdir(run):
            continue
        F = frames(run, keys=("points", "velocity", "alpha_s", "p_e", "phi_s"))
        ts, xf, pb, um = [], [], [], []
        for t, d in F:
            xs, hs = slide_height(d)
            ts.append(t); xf.append(front(xs, hs, 0.001))
            # excess pore pressure at the bed, x = 0.02 m (the sensor of Rondon et al.)
            P = d["points"]; x = np.round(P[:, 0], 9); z = P[:, 2]
            xv = np.unique(x)[np.argmin(abs(np.unique(x) - 0.02))]
            s = (x == xv) & (P[:, 1] == P[:, 1].min())
            zz, idx = np.unique(np.round(z[s], 9), return_index=True)
            pb.append(float(np.interp(0.5*(zz[1] - zz[0]), zz, d["p_e"][s][idx])))
            um.append(float(np.max(np.linalg.norm(d["velocity"], axis=1)[d["alpha_s"] > 0.5])) if np.any(d["alpha_s"] > 0.5) else 0.0)
        ts = np.array(ts); xf = np.array(xf)
        xs, hs = slide_height(F[-1][1])
        # end of the collapse: the front within 1 % of its final position
        ff = xf[-1] - L0
        t95 = float(ts[np.where(xf - L0 >= 0.95*ff)[0][0]]) if ff > 0 else 0.0
        rows[lab] = dict(t_end=round(float(ts[-1]), 2), runout_m=round(float(ff), 4), runout_L0=round(float(ff/L0), 3),
                         t95_s=round(t95, 2), H_end_over_H0=round(float(hs[0]/H0), 3), pe_base_min=round(float(min(pb)), 1),
                         pe_base_max=round(float(max(pb)), 1), umax=round(float(max(um)), 4), wall_s=wall(run))
        V = log_volume(run)
        if len(V):
            rows[lab]["volume_rel_change_max"] = float("%.2e" % np.max(abs(V[:, 1])))
        ax[0].plot(xs*100, hs*100, color=col, lw=1.6, label=lab)
        ax[0].plot([0, 6, 6], [H0*100, H0*100, 0], color=col, lw=0.8, ls=":")
        ax[1].plot(ts, (xf - L0)*100, color=col, lw=1.6, label=lab)
        ax[2].plot(ts, pb, color=col, lw=1.6, label=lab)
    if len(rows) == 2:
        k = list(rows)
        rows["runout dense/loose"] = round(rows[k[1]]["runout_m"]/rows[k[0]]["runout_m"], 3)
    ax[0].set_xlabel("x (cm)"); ax[0].set_ylabel("h (cm)"); ax[0].legend(fontsize=8); ax[0].set_title("deposit (end of run)", fontsize=10)
    ax[1].set_xlabel("t (s)"); ax[1].set_ylabel(r"$L - L_0$ (cm)"); ax[1].set_xscale("log"); ax[1].legend(fontsize=8)
    ax[2].set_xlabel("t (s)"); ax[2].set_ylabel(r"$p_e$ at the bed, x = 2 cm (Pa)"); ax[2].set_xscale("log"); ax[2].axhline(0, color=GREY, lw=0.8)
    fig.tight_layout(); fig.savefig(os.path.join(out, "14_rondon_loose_dense.png"))
    return rows


def robbe(runs, out):
    plt = R.style()
    h0 = 0.05
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
    rows = {}
    for k, (run, lab, col) in enumerate(zip(runs, ("dx 10 mm", "dx 5 mm"), (LBLUE, BLUE))):
        if run is None or not os.path.isdir(run):
            continue
        F = frames(run, keys=("points", "velocity", "phi", "alpha_s"))
        ts, Am, xc, xg = [], [], [], []
        snaps = []
        for t, d in F:
            x, e = water_surface(d)
            xs, hs = slide_height(d)
            # granular front at the water level: the furthest column with slide above z = 0 beyond the step
            xgf = front(xs, hs, 0.005)
            # ahead of the slide, before the far wall (x < 1.3 m: the reflection reaches there at about 1.2 s)
            w = (x > max(xgf, 0.0) + 0.02) & (x < 1.3) & np.isfinite(e)
            if np.any(w):
                i = np.argmax(e[w]); Am.append(float(e[w][i] - h0)); xc.append(float(x[w][i]))
            else:
                Am.append(np.nan); xc.append(np.nan)
            ts.append(t); xg.append(xgf); snaps.append((t, x, e))
        ts = np.array(ts); Am = np.array(Am); xc = np.array(xc)
        Am = np.where(ts < 1.15, Am, np.nan)
        im = int(np.nanargmax(Am))
        # crest speed after the wave has separated (t > 0.6 s): slope of the crest position
        m = (ts > 0.5) & (ts < 1.05) & np.isfinite(xc)
        c = float(np.polyfit(ts[m], xc[m], 1)[0]) if np.sum(m) > 3 else np.nan
        V = log_volume(run)
        rows[lab] = dict(A_max_cm=round(100*float(Am[im]), 2), t_A_max=round(float(ts[im]), 2), x_A_max=round(float(xc[im]), 3),
                         crest_speed=round(c, 3), granular_front_1s=round(float(np.interp(1.0, ts, xg)), 3), wall_s=wall(run),
                         volume_rel_change_max=float("%.2e" % np.max(abs(V[:, 1]))) if len(V) else None)
        ax[0].plot(ts, 100*Am, color=col, lw=1.6, label=lab)
        ax[1].plot(ts, xc, color=col, lw=1.6, label=lab)
        if k == len(runs) - 1 or runs[-1] is None:
            for (t, x, e), cc in zip([s for s in snaps if any(abs(s[0] - tt) < 0.011 for tt in (0.3, 0.6, 0.9, 1.2))], (GREY, GREEN, ORANGE, BLUE)):
                ax[2].plot(x, 100*(e - h0), color=cc, lw=1.4, label="t = %.1f s" % t)
    ax[0].axhline(8.25, color=GREY, ls="--", lw=1.0, label="exp.: 2.15 (V - 0.15)^0.58")
    ax[0].axvline(0.64, color=GREY, ls=":", lw=1.0)
    ax[0].set_xlabel("t (s)"); ax[0].set_ylabel("A (cm)"); ax[0].legend(fontsize=8); ax[0].set_title("wave amplitude", fontsize=10)
    tt = np.linspace(0.6, 1.5, 10)
    ax[1].plot(tt, 0.3 + 1.15*(tt - 0.6), "--", color=GREY, label="1.15 m/s (exp.)")
    ax[1].set_xlabel("t (s)"); ax[1].set_ylabel("crest x (m)"); ax[1].legend(fontsize=8)
    ax[2].set_xlabel("x (m)"); ax[2].set_ylabel(r"$\eta - h_0$ (cm)"); ax[2].legend(fontsize=8); ax[2].set_xlim(0, 1.6)
    fig.tight_layout(); fig.savefig(os.path.join(out, "15_robbe_saule_wave.png"))
    return rows


if __name__ == "__main__":
    cmd = sys.argv[1]
    a = [None if s == "-" else s for s in sys.argv[2:]]
    if cmd == "filmair":
        r = filmair(a[0], a[1], a[2], a[3], a[-1])
    elif cmd == "current":
        r = current(a[0], None, a[-1])
    elif cmd == "collapse":
        r = collapse(a[0], a[1], a[2], a[-1])
    elif cmd == "pelayer":
        r = pelayer([s for s in a[:-1] if s], a[-1])
    elif cmd == "rondon":
        r = rondon(a[0], a[1], a[-1])
    elif cmd == "robbe":
        r = robbe(a[:-1], a[-1])
    print(json.dumps(r, indent=1))
