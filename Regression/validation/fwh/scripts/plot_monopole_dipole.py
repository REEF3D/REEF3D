# Pressure fields of a harmonic point monopole and point dipole (linear acoustics, medium at rest),
# for comparison with the FW-H observer planes of validation case 15 (sphere, Re = 300, in water):
#
#   monopole:  p = A cos(w t - k r) / r
#   dipole:    p = A cos(theta) [ cos(w t - k r)/r^2 - k sin(w t - k r)/r ]     (axis y)
#
# The dipole has the hydrodynamic near field 1/r^2 for k r << 1 and the acoustic far field 1/r for
# k r >> 1; the change is at r = lambda/(2 pi). Figures (z = 0 plane):
#   1  snapshots near field (+-40 m) and acoustic field (+-40 km), scaled by r^2 or r
#   2  amplitude along the dipole axis against r, with the 1/r^2 and 1/r asymptotes
#   3  directivity (polar) of both sources
#   optional GIF animation of the acoustic field (--gif)
#
#   python plot_monopole_dipole.py [--f 0.137] [--c 1500] [--gif] [--noshow]     (or Run in Spyder)
# The figures are always written as PNG into the current directory; --noshow (or --save) only
# skips the windows.
import argparse, math, os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--f", type=float, default=0.137, help="frequency [Hz] (sphere Re 300: St U/D = 0.137)")
ap.add_argument("--c", type=float, default=1500.0, help="speed of sound [m/s]")
ap.add_argument("--gif", action="store_true", help="write dipole_far.gif (animation of the acoustic field)")
ap.add_argument("--noshow", "--save", action="store_true", help="do not show the figures (they are written anyway)")
args, _ = ap.parse_known_args()
if args.noshow:
    matplotlib.use("Agg")

w = 2*np.pi*args.f
k = w/args.c
lam = args.c/args.f
print(f"f = {args.f} Hz, c = {args.c} m/s: wavelength {lam:.0f} m, near/far change r = lambda/2pi = {1/k:.0f} m")


def monopole(x, y, t, A=1.0):
    r = np.hypot(x, y)
    return A*np.cos(w*t - k*r)/r


def dipole(x, y, t, A=1.0):
    # axis along y: cos(theta) = y/r
    r = np.hypot(x, y)
    return A*(y/r)*(np.cos(w*t - k*r)/r**2 - k*np.sin(w*t - k*r)/r)


def grid(L, n=400):
    # an even number of points: the grid does not pass through the singular centre r = 0
    s = np.linspace(-L, L, n)
    X, Y = np.meshgrid(s, s)
    return X, Y, np.hypot(X, Y)


# 1  snapshots: near field (+-40 m) scaled by r^2, acoustic field (+-40 km) scaled by r
fig1, ax = plt.subplots(2, 2, figsize=(11, 10))
for col, (L, scale, label) in enumerate(((40.0, 2, "near field, p r$^2$"), (40.0e3, 1, "acoustic field, p r"))):
    X, Y, R = grid(L)
    for row, (name, fun) in enumerate((("monopole", monopole), ("dipole (axis y)", dipole))):
        P = fun(X, Y, 0.0)*R**scale
        v = np.max(np.abs(P[R > 0.2*L]))     # colour range from outside the singular centre
        im = ax[row, col].pcolormesh(X/1e3 if L > 1e3 else X, Y/1e3 if L > 1e3 else Y, P,
                                     cmap="RdBu_r", vmin=-v, vmax=v, shading="auto", rasterized=True)
        ax[row, col].set_aspect("equal")
        unit = "km" if L > 1e3 else "m"
        ax[row, col].set_xlabel(f"x [{unit}]")
        ax[row, col].set_ylabel(f"y [{unit}]")
        ax[row, col].set_title(f"{name}: {label}")
        fig1.colorbar(im, ax=ax[row, col], shrink=0.8)
fig1.suptitle(f"t = 0, f = {args.f} Hz, c = {args.c} m/s, $\\lambda$ = {lam/1e3:.1f} km")
fig1.tight_layout()

# 2  dipole amplitude along its axis: 1/r^2 near, k/r far
r = np.logspace(0, 6, 400)
amp = np.sqrt(1/r**4 + k**2/r**2)
fig2, a2 = plt.subplots(figsize=(7, 5))
a2.loglog(r, amp, "k", lw=2, label="dipole |p| on the axis")
a2.loglog(r, 1/r**2, "--", label="near field 1/r$^2$")
a2.loglog(r, k/r, ":", label="acoustic field k/r")
a2.loglog(r, 1/r, color="gray", lw=1, label="monopole 1/r")
a2.axvline(1/k, color="r", lw=1)
a2.text(1/k, amp.max()*1e-2, " r = $\\lambda/2\\pi$", color="r")
a2.set_xlabel("r [m]")
a2.set_ylabel("|p| (A = 1)")
a2.set_ylim(1e-14, 2)
a2.grid(True, which="both", alpha=0.3)
a2.legend()
fig2.tight_layout()

# 3  directivity
th = np.linspace(0, 2*np.pi, 361)
fig3 = plt.figure(figsize=(9, 4.5))
for i, (name, d) in enumerate((("monopole", np.ones_like(th)), ("dipole (axis y)", np.abs(np.sin(th)))), 1):
    a3 = fig3.add_subplot(1, 2, i, projection="polar")
    a3.plot(th, d, lw=2)
    a3.set_title(name)
    a3.set_rticks([0.5, 1.0])
fig3.tight_layout()

# optional animation of the acoustic field of the dipole (one period)
if args.gif:
    from matplotlib.animation import FuncAnimation, PillowWriter
    X, Y, R = grid(40.0e3, 300)
    fig4, a4 = plt.subplots(figsize=(6, 5.5))
    P0 = dipole(X, Y, 0.0)*R
    v = np.max(np.abs(P0[R > 0.2*40.0e3]))
    im = a4.pcolormesh(X/1e3, Y/1e3, P0, cmap="RdBu_r", vmin=-v, vmax=v, shading="auto")
    a4.set_aspect("equal")
    a4.set_xlabel("x [km]")
    a4.set_ylabel("y [km]")
    fig4.colorbar(im, ax=a4, shrink=0.8, label="p r")
    nfr = 30

    def frame(i):
        t = i/(nfr*args.f)
        im.set_array((dipole(X, Y, t)*R).ravel())
        a4.set_title(f"dipole, acoustic field, t = {t:.2f} s")
        return (im,)

    FuncAnimation(fig4, frame, frames=nfr, blit=False).save("dipole_far.gif", writer=PillowWriter(fps=10))
    print("written dipole_far.gif")

fig1.savefig("monopole_dipole_fields.png", dpi=130)
fig2.savefig("dipole_near_far.png", dpi=130)
fig3.savefig("directivity.png", dpi=130)
print("written", ", ".join(os.path.abspath(f) for f in ("monopole_dipole_fields.png", "dipole_near_far.png", "directivity.png")))
if not args.noshow:
    plt.show()
