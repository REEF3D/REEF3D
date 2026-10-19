# P6c T7: the bed of a 180 deg bend against measurements: run 3 of Talmon & de Graaff (1989) in the curved flume of the
# Laboratory of Fluid Mechanics, TU Delft (data/talmon1989_run3.json): centre-line radius 4.10 m, width 0.50 m, straight
# reaches of 11.0 m upstream and 6.7 m downstream, Q 5.0 l/s, a0 0.051 m, slope 2.05 permille, Chezy 19.2 (ripples),
# D50 90 um, about 75 % of the transport in suspension.  SFLOW on a curvilinear grid with the total load of
# Engelund-Hansen (S 11 4) as bed load, the transverse bed slope effect (S 86, Talmon et al. 1995) and the spiral-flow
# correction (S 95); the first 2 m of the inflow reach keep their bed and feed the load (S 71, the sand supply of the
# flume 2 m below its entrance).  Compared: the relative depths a/a0 of the 48 cross-sections (0.32 m apart on the axis,
# section 45 the end of the bend) at 0.05 .. 0.45 m from the inner bank; a0 the mean depth of sections 1 - 5.
# Printed: rms of the model - measured a/a0 (all, the overshoot sections 6 - 24, the second half 25 - 45), the
# difference outer - inner (0.45 - 0.05 m) along the bend, the position and depth of the deepest outer-bank point.
# environment: KS (B 50), ESPIR (S 96, 1.0), MF (S 35, 100), DS (cell length, 0.025), STEPS (40000), XCTRL, NOSPIRAL, NOSLOPE,
# TOUT (the output interval, 10 s), TAVG (the averaging time, 300 s), TAG, NORUN, PNG
import sys, json, numpy as np
from p5lib import *
D = json.load(open(os.path.join(S5, "data", "talmon1989_run3.json")))
F, Rn = D["flume"], D["run"]
R, B, L1, L2 = F["Rc"], F["W"], F["L_in"], F["L_out"]
h0, Q, d50, Cz, slope = Rn["a0"], Rn["Q"], Rn["D50"], Rn["C"], Rn["i"]
ds = float(os.environ.get("DS", "0.025"))
steps = int(os.environ.get("STEPS", "60000"))
tout = float(os.environ.get("TOUT", "10"))           # the interval of the SFLOW output (P 182)
tavg = float(os.environ.get("TAVG", "300"))          # the bed (depth) averaged over the last TAVG s of flow
n1, n2, nb = int(round(L1/ds)), int(round(L2/ds)), int(round(np.pi*R/ds))
ny = int(round(B/ds)) | 1
nx = n1 + nb + n2
g, kappa, A, Bs, Cs = 9.81, 0.4, 9.0, 0.5, 0.3
espir = float(os.environ.get("ESPIR", "1.0"))
# the bed friction of SFLOW (A 218 1): Manning n = ks^(1/6)/20, so C = a0^(1/6)/n of the run gives ks (0.066 m: the
# whole resistance of the flume, ripples and side walls, on the bed; the run is then uniform at a0 0.051 m)
ks = float(os.environ.get("KS", (20.0*h0**(1.0/6.0)/Cz)**6))
# the bed shear of the sediment (S 16 2, Keulegan U/u* = 2.5 ln(12h/ks), ks = S 21 D50): the ks of C at a0
ks_sed = 12.0*h0/np.exp(Cz/np.sqrt(g)/2.5)

# nodes: s along the centre line, n to the left (towards the centre of the bend: positive n = inner bank)
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
Z = slope*(s_n[-1] - S)                            # the bed 0 at the outlet; F 65 2: the steady 1D flow from F 62

bc = "C 11 1\nC 12 21\nC 13 21\nC 14 2\nC 15 21\nC 16 3\n"
ctl = control(nx, ny, 1, nx*ds, ny*ds, 0.5, bc=bc)
ctrl = ("A 10 2\nA 210 3\nA 211 4\nA 218 1\nA 220 0\nA 240 1\nB 50 %g\nB 60 1\nW 10 %g\nF 60 %g\nF 62 %g\nF 65 2\nN 47 0.3\nM 10 1\nP 10 1\nP 181 -1\nP 182 %g\nW 22 -9.81\n"
        % (ks, Q, h0, h0, tout)
        + "S 10 1\nS 11 4\nS 16 2\nS 20 %g\nS 21 %g\nS 15 3\nS 17 1.0\nS 35 %s\nS 41 1\nS 43 6000\nS 44 1\nS 31 1\nS 100 0\nS 71 2.0\n"
        % (d50, ks_sed/d50, os.environ.get("MF", "100"))
        + "S 86 %d\nS 87 %g\nS 88 %g\nS 89 %g\nS 95 %d\nS 96 %g\n" % (0 if os.environ.get("NOSLOPE") else 1, A, Bs, Cs,
                                                                     0 if os.environ.get("NOSPIRAL") else 1, espir)
        + "G 50 2\n" + os.environ.get("XCTRL", "").replace(";", "\n") + "\n")
tag = os.environ.get("TAG", "")
d = os.path.join(S5, "runs/t7_talmon%s" % tag)
rc = 0
if not os.environ.get("NORUN"):
    rc = run(d, ctl, ctrl, grid=(X, Y, Z), steps=steps)
st = final_state(d)
print("rc=%d t=%.1f s  ks %.4f m  E_spir %g  MF %s" % (rc, st["simtime"], ks, espir, os.environ.get("MF", "100")),
      log_lines(d, "!!!")[:2])

WL = st["WL"].reshape(nx, ny)
# the ensemble average of the measurements (several soundings, the bed forms and the migrating bars averaged out): the
# depth of the SFLOW output (nodes) over the last TAVG s, mapped to (s, n) of the flume
import glob, re
sys.path.insert(0, os.path.join(S5, "../weber_junction"))
from weber_lib import read_vtp
def to_sn(x, y):
    s_, n_ = x.copy(), y.copy()
    up = (x <= L1) & (y < R)
    dn = (x <= L1) & (y >= R)
    s_[up], n_[up] = x[up], y[up]
    bb = x > L1
    r_, ph = np.hypot(x[bb] - L1, y[bb] - R), np.arctan2(y[bb] - R, x[bb] - L1)
    s_[bb], n_[bb] = L1 + R*(ph + 0.5*np.pi), R - r_
    s_[dn], n_[dn] = L1 + np.pi*R + (L1 - x[dn]), 2*R - y[dn]
    return s_, n_
vd = os.path.join(d, "REEF3D_SFLOW_VTP_FSF")
files = sorted(glob.glob(os.path.join(vd, "*.pvtp")) or glob.glob(os.path.join(vd, "*.vtp")))
acc, nav, tlast = None, 0, None
def vtp_time(fn):
    m = re.search(rb'Name="TimeValue"[^>]*>\s*([0-9.eE+-]+)', open(fn, "rb").read(4000))
    return float(m.group(1)) if m else None
tl = []
for fn in files:
    pcs = re.findall(r'Source="([^"]+)"', open(fn).read()) if fn.endswith(".pvtp") else [os.path.basename(fn)]
    V = [read_vtp(os.path.join(vd, pc)) for pc in pcs]
    t_ = vtp_time(os.path.join(vd, pcs[0]))
    tl.append((t_, V))
if tl:
    tend = max(t for t, _ in tl if t is not None)
    sel = [V for t, V in tl if t is not None and t >= tend - tavg]
    P0 = np.concatenate([v["Points"] for v in sel[0]])
    H = np.mean([np.concatenate([v["waterlevel"] for v in V]) for V in sel], axis=0)
    sP, nP = to_sn(P0[:, 0], P0[:, 1])
    from scipy.interpolate import LinearNDInterpolator
    fH = LinearNDInterpolator(np.c_[sP, nP], H)
    print("depth averaged over %d outputs, t %.0f .. %.0f s" % (len(sel), tend - tavg, tend))
else:
    fH = None
sc = 0.5*(s_n[:-1] + s_n[1:])
nc = 0.5*(n_n[:-1] + n_n[1:])
# the measured sections and lateral positions: s of section k on the axis, n = B/2 - (distance from the inner bank)
ks_ = np.arange(1, 49)
s_meas = L1 + np.pi*R + (ks_ - 45)*0.32
n_meas = 0.5*B - np.array(D["n_from_inner_m"])
def depth_at(s, n):
    i = np.clip(np.searchsorted(sc, s) - 1, 0, nx - 2); wi = np.clip((s - sc[i])/(sc[i+1] - sc[i]), 0, 1)
    j = np.clip(np.searchsorted(nc, n) - 1, 0, ny - 2); wj = np.clip((n - nc[j])/(nc[j+1] - nc[j]), 0, 1)
    return ((1-wi)*(1-wj)*WL[i, j] + wi*(1-wj)*WL[i+1, j] + (1-wi)*wj*WL[i, j+1] + wi*wj*WL[i+1, j+1])
A_mod = (np.array([[depth_at(s, n) for n in n_meas] for s in s_meas]) if fH is None
         else fH(np.array([[s, n] for s in s_meas for n in n_meas])).reshape(48, len(n_meas)))
a0m = A_mod[:5].mean()
A_mod = A_mod/a0m
A_mea = np.array(D["a_over_a0"])
dd = A_mod - A_mea
rms = lambda q: np.sqrt(np.mean(q**2))
print("a0 model %.4f m (measured %.3f)" % (a0m, h0))
print("a/a0: rms all %.3f   sections 6-24 %.3f   25-45 %.3f   (measured std over the bend %.3f)"
      % (rms(dd), rms(dd[5:24]), rms(dd[24:45]), A_mea[5:45].std()))
do_m, do_s = A_mea[:, -1] - A_mea[:, 0], A_mod[:, -1] - A_mod[:, 0]
print("outer - inner (0.45 - 0.05 m from the inner bank), mean of sections 25-45: measured %.2f model %.2f; sections 12-18: measured %.2f model %.2f"
      % (do_m[24:45].mean(), do_s[24:45].mean(), do_m[11:18].mean(), do_s[11:18].mean()))
km, kp = np.argmax(A_mea[:, -1]), np.argmax(A_mod[:, -1])
print("deepest at 0.45 m: measured section %d (%.2f a0)  model section %d (%.2f a0)" % (km+1, A_mea[km, -1], kp+1, A_mod[kp, -1]))
print("inner (0.05 m) minimum: measured section %d (%.2f a0)  model section %d (%.2f a0)"
      % (np.argmin(A_mea[:, 0])+1, A_mea[:, 0].min(), np.argmin(A_mod[:, 0])+1, A_mod[:, 0].min()))
np.savez(os.path.join(d, "t7.npz"), A_mod=A_mod, A_mea=A_mea, WL=WL, sc=sc, nc=nc)

png = os.environ.get("PNG")
if png:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.grid": True, "grid.color": "#e6e6e6"})
    fig = plt.figure(figsize=(12, 6.5))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1])
    lev = np.arange(0.2, 1.81, 0.1)
    for k, (Am, ttl) in enumerate(((A_mea, "measured (Talmon & de Graaff 1989, run 3)"), (A_mod, "SFLOW"))):
        ax = fig.add_subplot(gs[0, k])
        cs = ax.contourf(ks_, D["n_from_inner_m"], Am.T, levels=lev, cmap="RdYlBu", extend="both")
        ax.axvline(4.75, color="#6B6B6B", lw=0.8, ls="--"); ax.axvline(45, color="#6B6B6B", lw=0.8, ls="--")
        ax.set_title("a/a0, " + ttl, fontsize=9); ax.set_xlabel("cross-section (bend between the dashed lines)")
        ax.set_ylabel("distance from the inner bank (m)"); ax.grid(False)
    fig.colorbar(cs, ax=fig.axes[:2], shrink=0.8)
    ax = fig.add_subplot(gs[0, 2])
    ax.plot(ks_, do_m, "o", ms=3, mfc="none", mec="#6B6B6B", label="measured")
    ax.plot(ks_, do_s, color="#2B4C7E", lw=1.5, label="SFLOW")
    ax.axvline(4.75, color="#6B6B6B", lw=0.8, ls="--"); ax.axvline(45, color="#6B6B6B", lw=0.8, ls="--")
    ax.set_xlabel("cross-section"); ax.set_ylabel("(a outer - a inner)/a0  (0.45 and 0.05 m)")
    ax.legend(frameon=False)
    for k, cs_ in enumerate((10, 16, 22, 30, 40, 45)):
        ax = fig.add_subplot(2, 6, 7 + k)
        ax.plot(D["n_from_inner_m"], A_mea[cs_-1], "o", ms=3, mfc="none", mec="#6B6B6B")
        ax.plot(D["n_from_inner_m"], A_mod[cs_-1], color="#2B4C7E", lw=1.5)
        ang = (cs_ - 4.75)*0.32/R*180/np.pi
        ax.set_title("section %d (%.0f deg)" % (cs_, ang), fontsize=8); ax.set_ylim(0, 1.9)
        ax.set_xlabel("from the inner bank (m)")
        if k == 0: ax.set_ylabel("a/a0")
    fig.tight_layout()
    fig.savefig(png, dpi=150)
    print("wrote", png)
