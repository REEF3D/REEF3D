# Presentation figures of validation case 15 (sphere, Re = 300), in the style of the slides:
#   python sphere_slide_figs.py OUTDIR [t0 t1]      in the run directory (default window: second half)
# writes OUTDIR/forces.png, amplitude.png (+ amplitude.txt), directivity.png, timeseries.png.
# Needs numpy and matplotlib (font Helvetica Neue, else the matplotlib default).
import math, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = sys.argv[1] if len(sys.argv) > 1 else "slides"
os.makedirs(OUT, exist_ok=True)
DARK, LIGHT, BODY = "#0E2131", "#F5F3EE", "#3B4A57"
ORANGE, BLUE, GREY = "#D9662B", "#2C6E9E", "#8A96A0"
plt.rcParams.update({
    "font.family": "Helvetica Neue", "font.size": 15, "axes.titlesize": 16, "axes.labelsize": 16,
    "figure.facecolor": LIGHT, "axes.facecolor": LIGHT, "savefig.facecolor": LIGHT,
    "axes.edgecolor": BODY, "axes.labelcolor": DARK, "xtick.color": BODY, "ytick.color": BODY,
    "text.color": DARK, "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "legend.fontsize": 14, "lines.linewidth": 2.0})

rho, D = 1000.0, 1.0
c0, U, obs = 1500.0, 1.0, []
for l in open("ctrl.txt"):
    s = l.split()
    if len(s) >= 5 and s[0] == "U" and s[1] == "30":
        obs.append(np.array([float(x) for x in s[2:5]]))
q = 0.5*rho*U*U*math.pi*D*D/4.0


def table(fn, ncol):
    d = {}
    for l in open(fn):
        s = l.split()
        if len(s) < ncol:
            continue
        try:
            v = [float(x) for x in s[:ncol]]
        except ValueError:
            continue
        d[v[0]] = v[1:]
    t = np.array(sorted(d))
    return t, np.array([d[x] for x in t])


_cache = {}
def signal(n):
    if n not in _cache:
        d = {}
        for l in open(f"REEF3D_CFD_Acoustics/REEF3D-CFD-FWH-Observer-{n+1}.dat"):
            s = l.split()
            if len(s) != 2:
                continue
            try:
                d[float(s[0])] = float(s[1])
            except ValueError:
                pass
        t = np.array(sorted(d))
        _cache[n] = (t, np.array([d[x] for x in t]))
    return _cache[n]


tf, F6 = table("REEF3D_CFD_6DOF/REEF3D_6DOF_forces_0.dat", 13)
F6 = F6[:, 0:3]
tb, B = table("REEF3D_CFD_Acoustics/REEF3D-CFD-FWH-Box.dat", 13)
FB = -(B[:, 0:3] + B[:, 3:6] + np.gradient(B[:, 9:12], tb, axis=0))
dFB = np.gradient(FB, tb, axis=0)
d6 = np.gradient(F6, tf, axis=0)
tend = min(tf[-1], tb[-1])
T0 = float(sys.argv[2]) if len(sys.argv) > 2 else 0.5*tend
T1 = min(float(sys.argv[3]) if len(sys.argv) > 3 else tend, tend)
wb = (tb >= T0) & (tb <= T1)
mean = FB[wb].mean(axis=0)
ang = math.atan2(mean[2], mean[1])
# shedding frequency: peak of the lateral box force spectrum
tu = np.linspace(T0, T1, 4096)
lu = np.interp(tu, tb, FB[:, 1]*math.cos(ang) + FB[:, 2]*math.sin(ang))
spec = np.abs(np.fft.rfft((lu - lu.mean())*np.hanning(len(lu)), n=16*len(lu)))
fr = np.fft.rfftfreq(16*len(lu), tu[1]-tu[0])
sel = (fr > 0.02) & (fr < 1.0)
fst = fr[sel][np.argmax(spec[sel])]
k = 2*np.pi*fst/c0
print(f"window {T0:.1f}..{T1:.1f} s, St = {fst:.4f}")


def harmonic(t, p, f):
    w = 2*np.pi*f
    A = np.column_stack([np.ones_like(t), np.cos(w*t), np.sin(w*t)])
    c = np.linalg.lstsq(A, p, rcond=None)[0]
    return np.hypot(c[1], c[2]), np.arctan2(-c[2], c[1])


def curle(t, F, dF, x, tt):
    r = np.linalg.norm(x)
    rh = x/r
    tau = tt - r/c0
    Fi = np.column_stack([np.interp(tau, t, F[:, i]) for i in range(3)])
    dFi = np.column_stack([np.interp(tau, t, dF[:, i]) for i in range(3)])
    return -(Fi @ rh/r**2 + dFi @ rh/(c0*r))/(4*np.pi)


def win(n):
    t, p = signal(n)
    w = (t >= T0) & (t <= T1)
    return t[w], p[w]


# ---------- 1 forces
fig, ax = plt.subplots(2, 1, figsize=(13, 6.6), sharex=True)
for name, t, F, col in (("6DOF surface integral (direct forcing)", tf, F6, GREY), ("FW-H box momentum balance", tb, FB, BLUE)):
    m = t > 15
    ax[0].plot(t[m], F[m, 0]/q, lw=1.4, color=col, label=name)
    ax[1].plot(t[m], (F[m, 1]*math.cos(ang) + F[m, 2]*math.sin(ang))/q, lw=1.4, color=col)
ax[0].axhline(0.656, color=ORANGE, ls="--", lw=1.6, label="Johnson & Patel (1999)")
ax[1].axhline(0.069, color=ORANGE, ls="--", lw=1.6)
ax[0].set_ylim(0.5, 0.85)
ax[1].set_ylim(-0.02, 0.12)
for a in ax:
    a.axvspan(T0, T1, color="#E2E6E3", zorder=0)
ax[0].set_ylabel("C$_d$")
ax[1].set_ylabel("C$_l$ (shedding plane)")
ax[1].set_xlabel("t [s]")
ax[0].legend(loc="upper right", ncol=3, fontsize=13)
ax[1].text(T0 + 2, 0.105, f"averaging window {T0:.0f}–{T1:.0f} s", color=BODY, fontsize=13)
fig.tight_layout()
fig.savefig(f"{OUT}/forces.png", dpi=150)

# ---------- 2 amplitude against r (y axis)
iy, seen = [], set()
for n, x in enumerate(obs):
    if abs(x[0]) < 1e-6 and abs(x[2]) < 0.05 and x[1] > 1.0 and round(x[1], 3) not in seen:
        iy.append(n)
        seen.add(round(x[1], 3))
rows = []
for n in iy:
    x = obs[n]
    t, p = win(n)
    if len(t) < 50:
        continue
    t, p = t[::2], p[::2]
    rows.append((np.linalg.norm(x), harmonic(t, p, fst)[0], harmonic(t, curle(tb, FB, dFB, x, t), fst)[0],
                 harmonic(t, curle(tf, F6, d6, x, t), fst)[0]))
rows = np.array(sorted(rows))
Fy = harmonic(tb[wb], FB[wb, 1], fst)[0]
rr = np.logspace(0.3, 5.1, 300)
fig, (a, b) = plt.subplots(2, 1, figsize=(11, 8.2), sharex=True, gridspec_kw={"height_ratios": [3, 1.3]})
a.axvspan(rr[0], 1/k, color="#E8ECEE", zorder=0)
a.loglog(rr, Fy/(4*np.pi)*np.sqrt(1/rr**4 + k**2/rr**2), color=DARK, lw=1.6, label="compact dipole of the box force")
a.loglog(rows[:, 0], rows[:, 2], "s", mfc="none", mec=BLUE, mew=1.8, ms=11, label="Curle (box force)")
a.loglog(rows[:, 0], rows[:, 3], "^", mfc="none", mec=GREY, mew=1.5, ms=9, label="Curle (6DOF force)")
a.loglog(rows[:, 0], rows[:, 1], "o", color=ORANGE, ms=7, label="FW-H (permeable box)")
a.axvline(1/k, color=DARK, lw=1, ls=":")
a.text(1/k*1.12, 3e-2, f"r = λ/2π ≈ {1/k/1e3:.2f} km", color=DARK, fontsize=14)
a.text(2.6, 2e-9, "hydrodynamic near field  p ~ 1/r²", color=BODY, fontsize=14)
a.text(2.6e3, 3e-3, "acoustic far field\np ~ 1/r", color=BODY, fontsize=14)
a.set_ylabel(f"|p′| at f = {fst:.4f} Hz [Pa]")
a.legend(loc="upper right", bbox_to_anchor=(1.0, 0.80))
a.grid(True, which="major", alpha=0.25)
b.axhline(1, color=DARK, lw=1)
b.semilogx(rows[:, 0], rows[:, 1]/rows[:, 2], "o-", color=ORANGE, ms=7, lw=1.5)
b.set_ylim(0.9, 1.6)
for r, v in zip(rows[:, 0], rows[:, 1]/rows[:, 2]):
    if v > 1.6:
        b.annotate(f"{v:.1f}", (r, 1.58), ha="center", va="top", color=ORANGE, fontsize=13)
b.set_ylabel("FW-H / Curle")
b.set_xlabel("observer distance r on the y axis [m]")
b.grid(True, which="major", alpha=0.25)
fig.tight_layout()
fig.savefig(f"{OUT}/amplitude.png", dpi=150)
np.savetxt(f"{OUT}/amplitude.txt", rows, header="r FWH Curle_box Curle_6DOF")


# ---------- 3 directivity
def ring(R, tol):
    ids, seen = [], set()
    for n, x in enumerate(obs):
        a_ = round(math.degrees(math.atan2(x[1], x[0])), 1)
        if abs(np.hypot(x[0], x[1]) - R) < tol and abs(x[2]) < 0.05 and a_ not in seen:
            ids.append(n)
            seen.add(a_)
    return ids


fig = plt.figure(figsize=(12, 6))
for i, (lab, ids) in enumerate((("r = 20 km  (acoustic far field)", ring(2.0e4, 1.0)), ("r = 5 m  (2.5 m outside the box)", ring(5.0, 1e-3))), 1):
    th, af, ab = [], [], []
    for n in ids:
        t, p = win(n)
        t, p = t[::2], p[::2]
        th.append(math.atan2(obs[n][1], obs[n][0]))
        af.append(harmonic(t, p, fst)[0])
        ab.append(harmonic(t, curle(tb, FB, dFB, obs[n], t), fst)[0])
    o = np.argsort(th)
    th, af, ab = (np.append(np.array(v)[o], np.array(v)[o][0]) for v in (th, af, ab))
    ax3 = fig.add_subplot(1, 2, i, projection="polar")
    ax3.set_facecolor(LIGHT)
    s = ab.max()
    ax3.plot(th, ab/s, color=BLUE, lw=2.2, label="Curle (box force)")
    ax3.plot(th, af/s, "o", color=ORANGE, ms=7, label="FW-H")
    if i == 2:
        ax3.set_rscale("log")
        ax3.set_rlim(0.1, 30)
        ax3.set_rticks([0.3, 1, 3, 10])
        ax3.set_yticklabels(["0.3", "1", "3", "10"])
    else:
        ax3.set_rticks([0.5, 1.0])
    ax3.set_title(lab, pad=18)
    ax3.text(0, 0, "")
    ax3.grid(alpha=0.35)
    ax3.set_xticks(np.radians([0, 90, 180, 270]))
    ax3.set_xticklabels(["+x (flow)", "+y", "−x", "−y"])
    if i == 1:
        ax3.legend(loc="lower left", bbox_to_anchor=(-0.15, -0.12))
fig.tight_layout()
fig.savefig(f"{OUT}/directivity.png", dpi=150)

# ---------- 4 time series at 20 km and 100 m
fig, ax = plt.subplots(2, 1, figsize=(13, 6.6))
for a, ry in zip(ax, (2.0e4, 100.0)):
    n = min(iy, key=lambda j: abs(obs[j][1] - ry))
    t, p = signal(n)
    w = (t >= max(T0, T1 - 57.0)) & (t <= T1)
    t, p = t[w], p[w]
    cb = curle(tb, FB, dFB, obs[n], t)
    a.plot(t, (cb - cb.mean())*(1e9 if ry > 1e3 else 1e6), color=BLUE, lw=3.2, label="Curle (box force)")
    a.plot(t, (p - p.mean())*(1e9 if ry > 1e3 else 1e6), color=ORANGE, lw=1.5, label="FW-H")
    a.set_ylabel("p′ [nPa]" if ry > 1e3 else "p′ [µPa]")
    a.set_title(f"observer y = {ry/1e3:g} km  (delay r/c = {ry/c0:.1f} s)" if ry > 1e3 else f"observer y = {ry:g} m", loc="left")
ax[0].legend(loc="upper right", ncol=2)
ax[1].set_xlabel("observer time t [s]")
fig.tight_layout()
fig.savefig(f"{OUT}/timeseries.png", dpi=150)
print("figures done")
