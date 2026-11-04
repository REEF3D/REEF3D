# P7 T8: the lags of the bend morphology (sediment_exner_lag.cpp) in the LFM bend of T6 (SFLOW, curvilinear grid).
#   equiv    S 95 2 with S 97 0 (no adaptation length: I = I_e) against S 95 1 after STEPS steps with the bed moving:
#            the bed, the water level and the face velocities to round-off
#   profile  the bed held (S 35 0), the flow steady: the spiral-flow intensity I on the centre line (ST_spiral of the
#            bed output) with S 97 1 against the relaxation of the equilibrium I_e of the run with S 97 0, integrated
#            along the centre line with the scheme of the model (first-order upwind, backward: I_k = (I_e,k + c I_k-1)
#            / (1 + c), c = L_a/ds) and with the exact solution of I + L_a dI/ds = I_e; L_a = (m - 1) h/(2 kappa^2),
#            m = kappa C/sqrt(g) from the C and depth of the run (26.4, 0.2 m: 1.48 m).  Printed: the distance after the
#            bend entrance and exit to 63 % of the change (L_a of the model) and the rms of model - scheme
# (the lags across blocks and patches and the bed volume with S 64: t5_sed_mb.py with XCTRL "S 95 2;S 64 0.5")
# usage: python3 t8_bend_lags.py equiv|profile [steps]; environment DS (0.1), PNG, NORUN
import sys, glob, re, numpy as np
from p5lib import *
sys.path.insert(0, os.path.join(S5, "../weber_junction"))
from weber_lib import read_vtp

mode = sys.argv[1] if len(sys.argv) > 1 else "equiv"
steps = int(sys.argv[2]) if len(sys.argv) > 2 else (2000 if mode == "equiv" else 8000)
R, B, h0, Q, d50, Cz, slope = 4.25, 1.7, 0.2, 0.17, 0.00078, 26.4, 0.0018
L1, L2 = 11.0, 6.0
ds = float(os.environ.get("DS", "0.1"))
n1, n2, nb = int(round(L1/ds)), int(round(L2/ds)), int(round(np.pi*R/ds))
ny = int(round(B/ds)) | 1
nx = n1 + nb + n2
g, kappa = 9.81, 0.4
ks = 12.0*h0/np.exp(Cz/np.sqrt(g)/2.5)
ks_m = (20.0*h0**(1.0/6.0)/Cz)**6

s_n = np.concatenate([np.linspace(0, L1, n1+1), L1 + np.linspace(0, np.pi*R, nb+1)[1:], L1 + np.pi*R + np.linspace(0, L2, n2+1)[1:]])
n_n = np.linspace(-0.5*B, 0.5*B, ny+1)
S, N = np.meshgrid(s_n, n_n, indexing="ij")
X, Y = np.zeros_like(S), np.zeros_like(S)
m1, m2 = S <= L1, S >= L1 + np.pi*R
mb = ~(m1 | m2)
X[m1], Y[m1] = S[m1], N[m1]
phi = (S[mb] - L1)/R - 0.5*np.pi
r = R - N[mb]
X[mb], Y[mb] = L1 + r*np.cos(phi), R + r*np.sin(phi)
t = S[m2] - L1 - np.pi*R
X[m2], Y[m2] = L1 - t, 2*R - N[m2]
Z = slope*(s_n[-1] - S)

bc = "C 11 1\nC 12 21\nC 13 21\nC 14 2\nC 15 21\nC 16 3\n"
ctl = control(nx, ny, 1, nx*ds, ny*ds, 0.5, bc=bc)

def ctrl(extra, mf, out=False):
    return ("A 10 2\nA 210 3\nA 211 4\nA 218 1\nA 220 0\nA 240 1\nB 50 %g\nB 60 1\nW 10 %g\nF 60 %g\nF 62 %g\nF 65 2\nN 47 0.3\nM 10 1\nW 22 -9.81\n"
            % (ks_m, Q, h0, h0)
            + ("P 10 1\nP 76 1\nP 181 %d\nP 182 -1\n" % max(steps//4, 1) if out else "P 10 0\n")
            + "S 10 1\nS 11 2\nS 16 2\nS 20 %g\nS 21 %g\nS 15 3\nS 17 1.0\nS 35 %g\nS 41 1\nS 43 1000\nS 44 1\nS 31 1\nS 100 0\nS 71 2.0\n"
            % (d50, ks/d50, mf)
            + "S 86 1\nS 87 9\nS 88 0.5\nS 89 0.3\nS 96 0.4\nG 50 2\n" + extra.replace(";", "\n") + "\n")

def go(tag, c):
    d = os.path.join(S5, "runs/t8_%s" % tag)
    rc = 0
    if not os.environ.get("NORUN"):
        rc = run(d, ctl, c, grid=(X, Y, Z), steps=steps)
    st = final_state(d)
    print("%-10s rc=%d t=%.1f s" % (tag, rc, st["simtime"]), log_lines(d, "!!!")[:2])
    return d, st

def to_sn(x, y):
    s_, n_ = x.copy(), y.copy()
    up = (x <= L1) & (y < R)
    dn = (x <= L1) & (y >= R)
    bb = x > L1
    r_, ph = np.hypot(x[bb] - L1, y[bb] - R), np.arctan2(y[bb] - R, x[bb] - L1)
    s_[bb], n_[bb] = L1 + R*(ph + 0.5*np.pi), R - r_
    s_[dn], n_[dn] = L1 + np.pi*R + (L1 - x[dn]), 2*R - y[dn]
    return s_, n_

def centre(d, name):
    """the field of the last bed output on the centre line (mean of the node rows next to n = 0), at s of the cells"""
    vd = os.path.join(d, "REEF3D_SFLOW_VTP_BED")
    fs = sorted(glob.glob(os.path.join(vd, "*.pvtp")) or glob.glob(os.path.join(vd, "*.vtp")))
    fn = fs[-1]
    pcs = re.findall(r'Source="([^"]+)"', open(fn).read()) if fn.endswith(".pvtp") else [os.path.basename(fn)]
    V = [read_vtp(os.path.join(vd, pc)) for pc in pcs]
    P = np.concatenate([v["Points"] for v in V])
    f = np.concatenate([v[name] for v in V])
    sP, nP = to_sn(P[:, 0], P[:, 1])
    k = np.abs(nP) < 0.6*ds
    order = np.argsort(sP[k])
    sk, fk = sP[k][order], f[k][order]
    sc = 0.5*(s_n[:-1] + s_n[1:])
    # the two node rows at n = -ds/2, +ds/2 have the same s: average them, then to the cell centres
    su, inv = np.unique(np.round(sk, 6), return_inverse=True)
    fu = np.bincount(inv, fk)/np.bincount(inv)
    return sc, np.interp(sc, su, fu)

if mode == "equiv":
    _, a = go("equiv_s95_1", ctrl("S 95 1", 100))
    _, b = go("equiv_s95_2", ctrl("S 95 2;S 97 0", 100))
    for k in ("bed", "WL", "P", "Q"):
        print("S 95 2 / S 97 0 vs S 95 1: %-3s max diff %.2e (max |value| %.2e)" % (k, np.abs(a[k] - b[k]).max(), np.abs(a[k]).max()))

elif mode == "profile":
    m = kappa*Cz/np.sqrt(g)
    La = (m - 1.0)*h0/(2*kappa**2)
    d0, _ = go("prof_s97_0", ctrl("S 95 2;S 97 0", 0.0, out=True))
    d1, _ = go("prof_s97_1", ctrl("S 95 2;S 97 1", 0.0, out=True))
    sc, Ie = centre(d0, "ST_spiral")
    _, I1 = centre(d1, "ST_spiral")
    ds_c = np.diff(s_n)
    # the scheme of the model (backward, upwind) and the exact relaxation of I_e along s
    Is, Ix = np.zeros_like(Ie), np.zeros_like(Ie)
    for k in range(1, len(Ie)):
        c = La/ds_c[k]
        Is[k] = (Ie[k] + c*Is[k-1])/(1 + c)
        h_ = sc[k] - sc[k-1]
        e = np.exp(-h_/La)
        Ix[k] = Ix[k-1]*e + (1 - e)*0.5*(Ie[k] + Ie[k-1])     # exact for I_e linear between the cells
    Imax = np.abs(Ie).max()
    sb0, sb1 = L1, L1 + np.pi*R
    def reach63(I, s0, up):
        """distance after s0 at which I has made 63 % of its change from I(s0) to the level of I_e 3 m on"""
        k0 = np.searchsorted(sc, s0)
        k3 = np.searchsorted(sc, s0 + 3.0)
        a, b = I[k0], Ie[k3]
        tgt = a + (1 - np.exp(-1))*(b - a)
        kk = np.arange(k0, k3)
        hit = kk[(I[kk] - tgt)*np.sign(b - a) >= 0]
        return sc[hit[0]] - s0 if len(hit) else np.nan
    print("L_a nominal (C %.1f, h %.2f m) %.3f m; the spiral intensity |I_e| max %.4f m/s" % (Cz, h0, La, Imax))
    print("63 %% of the change after the bend entrance: model %.2f m, scheme %.2f m, exact %.2f m"
          % (reach63(I1, sb0, True), reach63(Is, sb0, True), reach63(Ix, sb0, True)))
    print("63 %% of the change after the bend exit:     model %.2f m, scheme %.2f m, exact %.2f m"
          % (reach63(I1, sb1, False), reach63(Is, sb1, False), reach63(Ix, sb1, False)))
    print("rms model - scheme %.3f, model - exact %.3f, I_e - exact %.3f (of |I_e| max)"
          % (np.sqrt(np.mean((I1 - Is)**2))/Imax, np.sqrt(np.mean((I1 - Ix)**2))/Imax, np.sqrt(np.mean((Ie - Ix)**2))/Imax))
    if os.environ.get("PNG"):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.spines.top": False,
                             "axes.spines.right": False, "axes.grid": True, "grid.color": "#e6e6e6"})
        fig, ax = plt.subplots(figsize=(8, 3.6))
        ax.plot(sc, Ie/Imax, color="#6B6B6B", lw=1.2, label="equilibrium $I_e = h|U|K_c$ (S 97 0)")
        ax.plot(sc, Ix/Imax, color="#4F8FCF", lw=2.4, alpha=0.6, label="exact relaxation, $L_a$ = %.2f m" % La)
        ax.plot(sc, Is/Imax, "--", color="#2B4C7E", lw=1.2, label="upwind scheme, $L_a$ = %.2f m" % La)
        ax.plot(sc, I1/Imax, color="#E07B39", lw=1.4, label="model (S 95 2, S 97 1)")
        for x_ in (sb0, sb1):
            ax.axvline(x_, color="#6B6B6B", lw=0.7, ls=":")
        ax.set_xlim(L1 - 3, sb1 + 5)
        ax.set_xlabel("s along the centre line (m); bend from %.1f to %.1f m" % (sb0, sb1))
        ax.set_ylabel("$I/|I_e|_{max}$")
        ax.legend(frameon=False, fontsize=8, loc="lower center")
        ax.set_title("Spiral-flow intensity with adaptation length, LFM bend (SFLOW, curvilinear, fixed bed)", fontsize=9)
        fig.tight_layout()
        fig.savefig(os.environ["PNG"], dpi=150)
        print("wrote", os.environ["PNG"])
