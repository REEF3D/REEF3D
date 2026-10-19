# Figures of a KP505 open-water run with FW-H for slides (light background):
#   python kp505_slide_figs.py OUTDIR [t0 t1]         in the run directory (default 0.4-0.8 s)
# kp505_bpf.png      BPF amplitude of FW-H and the CFD probes at each near point (bars) and the phase
# kp505_signal.png   p'(t) at 2R in the rotor plane, FW-H on the probe, and their spectra
# kp505_profiles.png time-mean axial velocity behind the rotor (from the LAGOON store)
import math, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

OUT = sys.argv[1] if len(sys.argv) > 1 else "slides"
T0 = float(sys.argv[2]) if len(sys.argv) > 2 else 0.4
T1 = float(sys.argv[3]) if len(sys.argv) > 3 else 0.8
os.makedirs(OUT, exist_ok=True)
DARK, LIGHT, BODY = "#0E2131", "#F5F3EE", "#3B4A57"
ORANGE, BLUE, GREY = "#D9662B", "#2C6E9E", "#8A96A0"
plt.rcParams.update({"font.family": "Helvetica Neue", "font.size": 15, "figure.facecolor": LIGHT, "axes.facecolor": LIGHT,
                     "savefig.facecolor": LIGHT, "axes.edgecolor": BODY, "axes.labelcolor": DARK, "xtick.color": BODY,
                     "ytick.color": BODY, "text.color": DARK, "axes.spines.top": False, "axes.spines.right": False,
                     "legend.frameon": False})

n = None; obs = []; npr = 0
for l in open("ctrl.txt"):
    s = l.split()
    if len(s) >= 3 and s[0] == "X" and s[1] == "211": n = abs(float(s[2]))/(2*math.pi)
    if len(s) >= 5 and s[0] == "U" and s[1] == "30": obs.append([float(v) for v in s[2:5]])
    if len(s) >= 5 and s[0] == "P" and s[1] == "64": npr += 1
fb = 5*n


def rows(fn):
    d = {}
    for l in open(fn):
        s = l.split()
        if len(s) == 2:
            try: d[float(s[0])] = float(s[1])
            except ValueError: pass
    t = np.array(sorted(d)); return t, np.array([d[x] for x in t])


def harm(t, p, f):
    A = np.column_stack([np.ones_like(t), np.cos(2*np.pi*f*t), np.sin(2*np.pi*f*t)])
    c = np.linalg.lstsq(A, p, rcond=None)[0]
    return math.hypot(c[1], c[2]), math.degrees(math.atan2(-c[2], c[1]))


labels, af, ap_, phf, php = [], [], [], [], []
for k in range(npr):
    x = obs[k]
    tf, pf = rows(f"REEF3D_CFD_Acoustics/REEF3D-CFD-FWH-Observer-{k+1}.dat")
    tp, pp = rows(f"REEF3D_CFD_PressureProbe/REEF3D-CFD-Probe-Pressure-{k+1}.dat")
    wf = (tf >= T0) & (tf <= T1); wp = (tp >= T0) & (tp <= T1)
    a1, p1 = harm(tf[wf], pf[wf], fb); a2, p2 = harm(tp[wp], pp[wp], fb)
    r = math.hypot(x[1], x[2])
    labels.append(f"x {x[0]:+.2f}\nr {r:.3f}")
    af.append(a1); ap_.append(a2); phf.append(p1); php.append(p2)
af, ap_ = np.array(af), np.array(ap_)

fig, ax = plt.subplots(figsize=(12, 5.2))
i = np.arange(len(af))
ax.bar(i - 0.2, ap_, 0.4, color=BLUE, label="CFD pressure probe")
ax.bar(i + 0.2, af, 0.4, color=ORANGE, label="FW-H (permeable box)")
for k in range(len(af)):
    dph = (phf[k] - php[k] + 180) % 360 - 180
    ax.text(k, max(af[k], ap_[k])*1.03 + 0.1, f"{100*(af[k]/ap_[k]-1):+.0f} %\n{dph:+.0f}°", ha="center", fontsize=12, color=BODY)
ax.set_xticks(i)
ax.set_xticklabels(labels, fontsize=12)
ax.set_ylabel(f"|p'| at BPF = {fb:.0f} Hz [Pa]")
ax.set_ylim(0, max(af.max(), ap_.max())*1.45)
ax.legend(loc="upper center", ncol=2)
fig.tight_layout()
fig.savefig(f"{OUT}/kp505_bpf.png", dpi=150)

# signal and spectrum at observer 1
tf, pf = rows("REEF3D_CFD_Acoustics/REEF3D-CFD-FWH-Observer-1.dat")
tp, pp = rows("REEF3D_CFD_PressureProbe/REEF3D-CFD-Probe-Pressure-1.dat")
fig, ax = plt.subplots(1, 2, figsize=(13, 4.8), gridspec_kw={"width_ratios": [1.6, 1]})
w = (tp >= T1 - 0.1) & (tp <= T1); wf = (tf >= T1 - 0.1) & (tf <= T1)
ax[0].plot(tp[w], pp[w] - pp[w].mean(), color=BLUE, lw=3.2, label="CFD probe")
ax[0].plot(tf[wf], pf[wf] - pf[wf].mean(), color=ORANGE, lw=1.6, label="FW-H")
ax[0].set_xlabel("t [s]  (one revolution)")
ax[0].set_ylabel("p' [Pa]")
ax[0].set_title(f"observer at 2R in the rotor plane {tuple(obs[0])}", loc="left", fontsize=14)
ax[0].legend(loc="upper right", ncol=2)
for (t, p), col, lab in (((tp, pp), BLUE, "CFD probe"), ((tf, pf), ORANGE, "FW-H")):
    ww = (t >= T0) & (t <= T1)
    tu = np.linspace(T0, T1, 4096); pu = np.interp(tu, t[ww], p[ww])
    sp = 2*np.abs(np.fft.rfft((pu - pu.mean())*np.hanning(len(pu))))/(0.5*len(pu))
    fr = np.fft.rfftfreq(len(pu), tu[1] - tu[0])
    ax[1].semilogy(fr, sp, color=col, lw=2.2 if col == BLUE else 1.4, label=lab)
for m in (1, 2, 3):
    ax[1].axvline(m*fb, color=GREY, lw=0.8, ls=":")
ax[1].set_xlim(0, 4*fb)
ax[1].set_ylim(1e-3, None)
ax[1].set_xlabel("f [Hz]  (BPF and harmonics dotted)")
ax[1].set_ylabel("amplitude [Pa]")
fig.tight_layout()
fig.savefig(f"{OUT}/kp505_signal.png", dpi=150)

# axial velocity profiles from the store
if os.path.isdir("REEF3D_CFD.lagoon"):
    from lagoon_read import lagoon
    S = lagoon("REEF3D_CFD.lagoon"); x, y, z = S.coords(); Y, Z = np.meshgrid(y, z); R = np.hypot(Y, Z)
    fr = [k for k, t in enumerate(S.times) if T0 + 0.2 <= t <= T1]
    rb = np.arange(0.0, 0.25, 0.005)
    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    for xv, col in ((-0.05, GREY), (0.05, ORANGE), (0.15, BLUE), (0.3, DARK)):
        ii = int(np.argmin(abs(x - xv))); acc = 0
        for k in fr:
            acc = acc + S.volume("velocity", k)[:, :, ii, 0]
        ux = acc/len(fr)
        prof = [ux[(R >= a) & (R < a + 0.005)].mean() for a in rb]
        ax.plot(np.array(prof), rb + 0.0025, color=col, lw=2.2, label=f"x = {x[ii]:+.2f} m")
    ax.axhline(0.125, color=GREY, lw=0.8, ls="--")
    ax.text(0.05, 0.128, "R", color=BODY)
    ax.set_xlabel("time-mean axial velocity u [m/s]  (U = 1.5)")
    ax.set_ylabel("r [m]")
    ax.legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(f"{OUT}/kp505_profiles.png", dpi=150)
print("written", OUT)
