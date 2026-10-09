# Figures of validation case 15 (sphere, Re = 300), in the run directory:
#   python sphere_plots.py [t0] [t1] [--noshow]     (default window: the second half of the output)
# The figures are always written as sphere_fig1..4.png into the run directory; --noshow (or the
# old --save) only skips the windows, e.g. on a machine without a display.
#   1  forces: Cd and the lateral Cl (6DOF surface integral and FW-H box momentum balance)
#   2  amplitude at the shedding frequency along the y axis against r: FW-H, Curle (box force and
#      6DOF force) and the dipole of the box force |F_y| sqrt(1/r^4 + k^2/r^2)/(4 pi)
#   3  directivity on the rings r = 5 m and r = 20 km: FW-H and Curle (box force)
#   4  time series at y = 100 m and y = 20 km: FW-H and Curle (box force)
# Needs numpy and matplotlib; reads ctrl.txt, REEF3D_CFD_6DOF/, REEF3D_CFD_Acoustics/.
import math, os, sys
import numpy as np
import matplotlib

args = [a for a in sys.argv[1:] if not a.startswith("--")]
noshow = "--noshow" in sys.argv or "--save" in sys.argv
if noshow:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt

rho, D = 1000.0, 1.0
c0, U, nu, obs = 1500.0, None, None, []
for l in open("ctrl.txt"):
    s = l.split()
    if len(s) < 3:
        continue
    if s[0] == "U" and s[1] == "11": c0 = float(s[2])
    if s[0] == "U" and s[1] == "50": U = float(s[2])
    if s[0] == "W" and s[1] == "2": nu = float(s[2])
    if s[0] == "U" and s[1] == "30": obs.append(np.array([float(x) for x in s[2:5]]))
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


def signal(n):
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
    return t, np.array([d[x] for x in t])


tf, F6 = table("REEF3D_CFD_6DOF/REEF3D_6DOF_forces_0.dat", 13)
F6 = F6[:, 0:3]
tb, B = table("REEF3D_CFD_Acoustics/REEF3D-CFD-FWH-Box.dat", 13)
dM = np.gradient(B[:, 9:12], tb, axis=0)
FB = -(B[:, 0:3] + B[:, 3:6] + dM)

tend = min(tf[-1], tb[-1])
T0 = float(args[0]) if len(args) > 0 else 0.5*tend
T1 = min(float(args[1]) if len(args) > 1 else tend, tend)
if np.count_nonzero((tf >= T0) & (tf <= T1)) < 50:
    sys.exit(f"window {T0}..{T1} s is outside the output (t <= {tend:.2f} s)")
print(f"Re = {U*D/nu:.0f}, window {T0:.1f}..{T1:.1f} s")


def harmonic(t, p, f):
    # least-squares mean, amplitude and phase of the component at f
    w = 2*np.pi*f
    A = np.column_stack([np.ones_like(t), np.cos(w*t), np.sin(w*t)])
    c = np.linalg.lstsq(A, p, rcond=None)[0]
    return np.hypot(c[1], c[2]), np.arctan2(-c[2], c[1])


# shedding frequency from the lateral box force (uniform resampling, FFT with zero padding)
m = (tb >= T0) & (tb <= T1)
mean = FB[m].mean(axis=0)
ang = math.atan2(mean[2], mean[1])
lat = FB[:, 1]*math.cos(ang) + FB[:, 2]*math.sin(ang)
tu = np.linspace(T0, T1, 4096)
lu = np.interp(tu, tb, lat)
spec = np.abs(np.fft.rfft((lu - lu.mean())*np.hanning(len(lu)), n=16*len(lu)))
fr = np.fft.rfftfreq(16*len(lu), tu[1]-tu[0])
sel = (fr > 0.02) & (fr < 1.0)
fst = fr[sel][np.argmax(spec[sel])]
k = 2*np.pi*fst/c0
for name, t, F in (("6DOF surface integral", tf, F6), ("box momentum balance", tb, FB)):
    w = (t >= T0) & (t <= T1)
    mF = F[w].mean(axis=0)
    l = F[w, 1]*math.cos(ang) + F[w, 2]*math.sin(ang)
    a, _ = harmonic(t[w], l, fst)
    print(f"{name:>22}: Cd = {mF[0]/q:.4f}, mean Cl = {np.hypot(mF[1], mF[2])/q:.4f}, Cl' = {a/q:.4f}")
print(f"St = {fst*D/U:.4f}   (Johnson & Patel 1999: Cd 0.656, Cl 0.069, St 0.137)")


def curle(t, F, x, tt):
    r = np.linalg.norm(x)
    rh = x/r
    tau = tt - r/c0
    Fi = np.column_stack([np.interp(tau, t, F[:, i]) for i in range(3)])
    dF = np.column_stack([np.interp(tau, t, np.gradient(F[:, i], t)) for i in range(3)])
    return -(Fi @ rh/r**2 + dF @ rh/(c0*r))/(4*np.pi)


def observer(n, x):
    t, p = signal(n)
    w = (t >= T0) & (t <= T1)
    return t[w], p[w]


# observer groups
iy, seen = [], set()
for n, x in enumerate(obs):
    # the observers on the y axis; the same points can come from the line, the ring and the planes
    if abs(x[0]) < 1e-6 and abs(x[2]) < 0.05 and x[1] > 1.0 and round(x[1], 3) not in seen:
        iy.append(n)
        seen.add(round(x[1], 3))
def ring(R, tol):
    ids, seen = [], set()
    for n, x in enumerate(obs):
        a = round(math.degrees(math.atan2(x[1], x[0])), 1)
        if abs(np.hypot(x[0], x[1]) - R) < tol and abs(x[2]) < 0.05 and a not in seen:
            ids.append(n)
            seen.add(a)
    return ids
ring5 = ring(5.0, 1e-3)
ring20 = ring(2.0e4, 1.0)

# 1 forces
fig1, ax1 = plt.subplots(2, 1, figsize=(9, 6), sharex=True)
for name, t, F in (("6DOF", tf, F6), ("box", tb, FB)):
    ax1[0].plot(t, F[:, 0]/q, lw=0.8, label=name)
    ax1[1].plot(t, (F[:, 1]*math.cos(ang) + F[:, 2]*math.sin(ang))/q, lw=0.8, label=name)
# axis limits without the start-up (impulsive start, kick)
ts = 0.1*tend
vals = [np.concatenate([F[t > ts, 0]/q for t, F in ((tf, F6), (tb, FB))]),
        np.concatenate([(F[t > ts, 1]*math.cos(ang) + F[t > ts, 2]*math.sin(ang))/q for t, F in ((tf, F6), (tb, FB))])]
for a, lab, v in zip(ax1, ("Cd", "Cl (shedding plane)"), vals):
    lo, hi = np.percentile(v, 0.5), np.percentile(v, 99.5)
    a.set_ylim(lo - 0.2*(hi-lo), hi + 0.2*(hi-lo))
    a.axvspan(T0, T1, color="0.9")
    a.set_ylabel(lab)
    a.legend()
ax1[0].axhline(0.656, color="k", ls=":", lw=1)
ax1[1].set_xlabel("t [s]")
fig1.tight_layout()

# 2 amplitude against r on the y axis
rows = []
for n in iy:
    x = obs[n]
    t, p = observer(n, x)
    if len(t) < 20:
        continue
    r = np.linalg.norm(x)
    rows.append((r, harmonic(t, p, fst)[0], harmonic(t, curle(tb, FB, x, t), fst)[0],
                 harmonic(t, curle(tf, F6, x, t), fst)[0]))
rows = np.array(sorted(rows))
wb = (tb >= T0) & (tb <= T1)
Fy = harmonic(tb[wb], FB[wb, 1], fst)[0]
rr = np.logspace(0, 5.3, 300)
fig2, a2 = plt.subplots(figsize=(7.5, 5.5))
a2.loglog(rr, Fy/(4*np.pi)*np.sqrt(1/rr**4 + k**2/rr**2), "k", lw=1.5, label="dipole of the box force")
a2.loglog(rows[:, 0], rows[:, 2], "s", mfc="none", ms=9, label="Curle (box force)")
a2.loglog(rows[:, 0], rows[:, 3], "^", mfc="none", ms=8, label="Curle (6DOF force)")
a2.loglog(rows[:, 0], rows[:, 1], "o", ms=5, label="FW-H")
a2.axvline(1/k, color="r", lw=1)
a2.text(1/k, rows[:, 1].max()*0.3, " $\\lambda/2\\pi$", color="r")
a2.set_xlabel("r on the y axis [m]")
a2.set_ylabel(f"|p'| at f = {fst:.3f} Hz [Pa]")
a2.grid(True, which="both", alpha=0.3)
a2.legend()
fig2.tight_layout()
print("\n   r [m]     FW-H      Curle box  Curle 6DOF  FW-H/box")
for r, f, b, s6 in rows:
    print(f"{r:9.3g} {f:10.3e} {b:10.3e} {s6:10.3e}  {f/b:7.3f}")

# 3 directivity
rings = [(lab, ids) for lab, ids in (("r = 5 m", ring5), ("r = 20 km", ring20)) if ids]
fig3 = plt.figure(figsize=(5*len(rings), 4.8))
for i, (lab, ids) in enumerate(rings, 1):
    th, af, ab = [], [], []
    for n in ids:
        x = obs[n]
        t, p = observer(n, x)
        if len(t) < 20:
            continue
        th.append(math.atan2(x[1], x[0]))
        af.append(harmonic(t, p, fst)[0])
        ab.append(harmonic(t, curle(tb, FB, x, t), fst)[0])
    o = np.argsort(th)
    th, af, ab = (np.append(np.array(v)[o], np.array(v)[o][0]) for v in (th, af, ab))
    a3 = fig3.add_subplot(1, len(rings), i, projection="polar")
    a3.plot(th, ab/ab.max(), "k-", lw=1.5, label="Curle (box)")
    a3.plot(th, af/ab.max(), "o", ms=4, label="FW-H")
    a3.set_title(lab)
    a3.legend(loc="lower left", fontsize=8)
fig3.tight_layout()

# 4 time series
fig4, ax4 = plt.subplots(2, 1, figsize=(9, 6))
for a, ry in zip(ax4, (100.0, 2.0e4)):
    n = min(iy, key=lambda j: abs(obs[j][1] - ry))
    t, p = observer(n, obs[n])
    if len(t) < 20:
        continue
    a.plot(t, p - p.mean(), lw=1.2, label="FW-H")
    cb = curle(tb, FB, obs[n], t)
    a.plot(t, cb - cb.mean(), "--", lw=1.2, label="Curle (box)")
    a.set_title(f"observer y = {obs[n][1]:g} m")
    a.set_ylabel("p' [Pa]")
    a.legend()
ax4[1].set_xlabel("t [s]")
fig4.tight_layout()

for i, f in enumerate((fig1, fig2, fig3, fig4), 1):
    f.savefig(f"sphere_fig{i}.png", dpi=130)
print("written", ", ".join(os.path.abspath(f"sphere_fig{i}.png") for i in range(1, 5)))
if not noshow:
    plt.show()
