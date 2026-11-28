#!/usr/bin/env python3
# Architect: Hans Bihs
"""Evaluation of case 35 (dilatancy of the slide layer, W 164): a uniform submerged layer on an incline.

  dilatancy_eval.py <run dir> <out dir>

<run dir> holds the runs loose, dense, coulomb (SFLOW) and dense_nhflow (NHFLOW, optional). The plateau of the
layer (|h - h0| < 2 % h0) is a 0D problem: its speed u, solid fraction phi and tan(psi) are compared with the
ODE of the same closure (ode() below, small explicit / backward-Euler steps):

  du/dt   = (1 - r) g sin(theta) - max(mu + tan(psi), 0) max(sigma_0 - p_e, 0.02 sigma_0)/(rho_s h)
  dphi/dt = - phi gamma tan(psi),  gamma = c u/h,  p_e = - eta_f gamma tan(psi) h^2/(2 k)
  tan(psi) = K3 (phi - phi_c/(1 + J^1/2)),  J = eta_f gamma/max(sigma_0 - p_e, 0.02 sigma_0)

Reads the SFLOW/NHFLOW free-surface VTP files (P 18 2). Needs numpy and matplotlib.
"""
import sys, os, math, json
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import slide_eval as se

se.KEYS = se.KEYS + ("slide_phi", "slide_tanpsi")

# the parameters of the case (ctrl.txt)
P = dict(h=0.005, theta=25.0, delta=22.0, rhos=1870.0, rhow=998.2, d=1.6e-4, K3=4.09, phic=0.585, cg=3.0,
         etaf=998.2*1.004e-6)
PHI0 = dict(loose=0.55, dense=0.598, coulomb=0.598, dense_nhflow=0.598)


def ode(phi0, T, dil=True, dt=2.0e-4, h=0.005, theta=25.0, delta=22.0, rhos=1870.0, rhow=998.2, d=1.6e-4,
        K3=4.09, phic=0.585, cg=3.0, etaf=998.2*1.004e-6):
    g = 9.81
    th = math.radians(theta)
    gs = (1.0 - rhow/rhos)*g*math.sin(th)
    gn = g*math.cos(th)
    mu = math.tan(math.radians(delta))
    s0 = (rhos - rhow)*gn*h

    def perm(ph):
        ph = min(max(ph, 0.3), 0.7)
        return d*d*(1.0 - ph)**3/(180.0*ph*ph)

    def tps(u, ph):
        C1 = etaf*cg*h/(2.0*perm(ph))
        gam = cg*u/h
        if not dil:
            return 0.0, 0.0
        if gam <= 0.0:
            return min(max(K3*(ph - phic), -0.6), 0.6), C1
        C = C1*u
        f = lambda t: t - K3*(ph - phic/(1.0 + math.sqrt(etaf*gam/max(s0 + C*t, 0.02*s0))))
        a, b = -0.6, 0.6
        if f(a) >= 0.0:
            return a, C1
        if f(b) <= 0.0:
            return b, C1
        for _ in range(60):
            m = 0.5*(a + b)
            if f(m) > 0.0:
                b = m
            else:
                a = m
        return 0.5*(a + b), C1

    u, ph, t = 0.0, phi0, 0.0
    n = int(round(T/dt))
    every = max(1, int(round(0.02/dt)))
    out = []
    for i in range(n + 1):
        tp, C1 = tps(u, ph)
        if i % every == 0:
            out.append((t, u, ph, tp))
        if i == n:
            break
        m = max(mu + tp, 0.0)
        if tp > 0.0:
            # the suction is linear in the speed: backward Euler
            u = max(u + dt*(gs - m*s0/(rhos*h)), 0.0)/(1.0 + dt*m*C1*tp/(rhos*h))
        else:
            sp = max(s0 + C1*u*tp, 0.02*s0)
            u = max(u + dt*(gs - m*sp/(rhos*h)), 0.0)
        ph = min(max(ph*math.exp(-dt*cg*u/h*tp), 0.3), 0.75)
        t += dt
    return np.array(out)


def plateau(run, h0=0.005):
    """time series of the median speed, phi, tan(psi) of the plateau of the layer"""
    F = se.frames(run)
    rows = []
    for t, d in F:
        x, hh = se.line(d, "slide_h")
        _, u = se.line(d, "slide_velocity", comp=0)
        sel = abs(hh - h0) < 0.02*h0
        if sel.sum() < 20:
            continue
        r = [t, float(np.median(u[sel])), np.nan, np.nan, float(x[sel].min()), float(x[sel].max())]
        if "slide_phi" in d:
            r[2] = float(np.median(se.line(d, "slide_phi")[1][sel]))
            r[3] = float(np.median(se.line(d, "slide_tanpsi")[1][sel]))
        rows.append(r)
    return np.array(rows)


def trigger(t, u, uc=0.05):
    w = np.where(u > uc)[0]
    return float(t[w[0]]) if len(w) else float("nan")


def main(rundir, out):
    plt = se.style()
    os.makedirs(out, exist_ok=True)
    T = 15.0
    col = dict(loose=se.ORANGE, dense=se.BLUE, coulomb=se.GREY, dense_nhflow=se.GREEN)
    res = {}
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
    for c in ("loose", "dense", "coulomb", "dense_nhflow"):
        run = os.path.join(rundir, c)
        if not os.path.isdir(run):
            continue
        O = ode(PHI0[c], T, dil=(c != "coulomb"))
        M = plateau(run)
        uo = np.interp(M[:, 0], O[:, 0], O[:, 1])
        po = np.interp(M[:, 0], O[:, 0], O[:, 2])
        big = uo > 0.01
        r = dict(trigger_model=trigger(M[:, 0], M[:, 1]), trigger_ode=trigger(O[:, 0], O[:, 1]),
                 u_end_model=float(M[-1, 1]), u_end_ode=float(np.interp(M[-1, 0], O[:, 0], O[:, 1])),
                 t_end=float(M[-1, 0]),
                 u_rel_err_max=float(np.max(abs(M[big, 1] - uo[big])/uo[big])) if big.any() else float("nan"))
        if c != "coulomb":
            r["phi_end_model"] = float(M[-1, 2])
            r["phi_end_ode"] = float(np.interp(M[-1, 0], O[:, 0], O[:, 2]))
            r["phi_abs_err_max"] = float(np.nanmax(abs(M[:, 2] - po)))
        res[c] = r
        lab = c.replace("_", " ")
        if c != "dense_nhflow":
            ax[0].semilogy(O[:, 0], np.maximum(O[:, 1], 1e-7), color=col[c], lw=1.2, label=f"{lab}: ODE")
        ax[0].semilogy(M[:, 0], np.maximum(M[:, 1], 1e-7), "o" if c != "dense_nhflow" else "x", ms=2.5,
                       color=col[c], label=f"{lab}: {se.model(run)}")
        if c != "coulomb":
            if c != "dense_nhflow":
                ax[1].plot(O[:, 0], O[:, 2], color=col[c], lw=1.2)
                ax[2].plot(O[:, 0], O[:, 3], color=col[c], lw=1.2)
            ax[1].plot(M[:, 0], M[:, 2], "o" if c != "dense_nhflow" else "x", ms=2.5, color=col[c])
            ax[2].plot(M[:, 0], M[:, 3], "o" if c != "dense_nhflow" else "x", ms=2.5, color=col[c])
    ax[0].set(xlabel="t [s]", ylabel="u of the plateau [m/s]", ylim=(1e-6, 10))
    ax[0].legend(fontsize=7, ncol=2)
    ax[1].axhline(P["phic"], color=se.GREY, lw=0.6, ls=":")
    ax[1].set(xlabel="t [s]", ylabel=r"solid fraction $\phi$")
    ax[2].axhline(0.0, color=se.GREY, lw=0.6)
    ax[2].set(xlabel="t [s]", ylabel=r"tan $\psi$")
    fig.suptitle(r"35: submerged layer h = 5 mm on 25$^\circ$, Coulomb 22$^\circ$, d = 0.16 mm: loose $\phi_0$ = 0.55, dense $\phi_0$ = 0.598",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "35_dilatancy.png"))
    json.dump(res, open(os.path.join(out, "35_dilatancy.json"), "w"), indent=1)
    return res


if __name__ == "__main__":
    print(json.dumps(main(sys.argv[1], sys.argv[2]), indent=1))
