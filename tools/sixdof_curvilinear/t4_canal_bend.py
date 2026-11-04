# P5 gate: a floating box drifting with the current through a canal with a 90 deg bend (NHFLOW, X 10 1), on three
# blocks and on one block of the same nodes.
#   the canal (W 3 m, h 0.5 m, Q 0.6 m3/s, U 0.4 m/s): the inlet reach along +y (x 0 .. W, y -Lin .. 0), the bend to
#   the right about (W + R0, 0) (inner radius R0 7.5 m, the centreline radius 3 W), the outlet reach along +x
#   (y R0 .. R0 + W); i along the canal, j across from the inner bank (right-handed frames, the same in every block)
#   three: block 0 the outlet reach (native outflow, F 62), block 1 the inlet reach (inflow G 54 1 1), block 2 the bend
#   one:   one block of all nodes (native inflow B 60 1 and outflow)
# The box (0.5 x 0.5 x 0.3 m, density 500: draft 0.15 m) starts on the axis 3 m below the inlet; surge, sway, heave
# and yaw free; non-hydrostatic (with A 520 0 the heave of the box grew in the bend).  The method has no contact
# model: a body that reaches a bank or the outflow breaks the run (a box of 0.8 x 0.5 m in a canal of 2 m, centreline
# radius 2.5 W, kept its heading through the bend, lay across the flow and hit the outer bank), so the box is square,
# the canal wide and the path is evaluated up to 3 m before the outlet.  Printed: the path in canal coordinates (s
# along the axis from the start of the bend, n from the axis, positive towards the outer bank), the speed, the
# heading, three vs one.
import sys, os, glob, numpy as np
T6 = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(T6, "../nhflow_multiblock"))
from p4lib import run4

Wc, h0, U0 = 3.0, 0.5, 0.4
Q = U0*Wc*h0
dx = 0.125
n = int(round(Wc/dx))
R0 = 7.5
Rc = R0 + 0.5*Wc
Lin, Lout = 10.0, 12.0
ni_in, ni_out = int(round(Lin/dx)), int(round(Lout/dx))
ni_b = int(round(0.5*np.pi*Rc/dx))
nz = 5
tend = float(os.environ.get("TEND", "100"))
np_ = int(os.environ.get("NP", "1"))

def reach_in():
    I, J = np.meshgrid(np.arange(ni_in+1), np.arange(n+1), indexing="ij")
    return Wc - J*dx, -Lin + I*dx                         # i +y, j -x (j 0 at x = W, the inner bank)

def bend():
    I, J = np.meshgrid(np.arange(ni_b+1), np.arange(n+1), indexing="ij")
    th = np.pi - 0.5*np.pi*I/ni_b
    r = R0 + J*dx
    return Wc + R0 + r*np.cos(th), r*np.sin(th)

def reach_out():
    I, J = np.meshgrid(np.arange(ni_out+1), np.arange(n+1), indexing="ij")
    return Wc + R0 + I*dx, R0 + J*dx                     # i +x, j +y (j 0 at y = R0, the inner bank)

def blk(XY):
    X, Y = XY
    return X, Y, np.zeros_like(X)

B_in, B_b, B_out = blk(reach_in()), blk(bend()), blk(reach_out())
one = tuple(np.concatenate([B_in[k], B_b[k][1:], B_out[k][1:]], axis=0) for k in range(3))

def canal(x, y):
    """s along the axis from the start of the bend (negative in the inlet reach), n from the axis towards the outer
    bank, the direction of the axis"""
    x, y = np.asarray(x, float), np.asarray(y, float)
    s, nn, ang = np.empty_like(x), np.empty_like(x), np.empty_like(x)
    cx = Wc + R0
    for k in range(len(x)):
        if y[k] <= 0.0:
            s[k], nn[k], ang[k] = y[k], 0.5*Wc - x[k], 0.5*np.pi
        elif x[k] >= cx:
            s[k], nn[k], ang[k] = 0.5*np.pi*Rc + (x[k] - cx), y[k] - Rc, 0.0
        else:
            r, th = np.hypot(x[k]-cx, y[k]), np.arctan2(y[k], x[k]-cx)
            s[k], nn[k], ang[k] = (np.pi - th)*Rc, r - Rc, th - 0.5*np.pi
    return s, nn, ang

bxl, byl, bzl, rho = 0.5, 0.5, 0.3, 500.0
cx0, cy0 = 0.5*Wc, -Lin + 3.0
zb0 = h0 - rho/1000.0*bzl
body = ("X 10 1\nX 11 1 1 1 0 0 1\nX 21 %g\nX 110 %g %g %g %g %g %g\n"
        % (rho, cx0-0.5*bxl, cx0+0.5*bxl, cy0-0.5*byl, cy0+0.5*byl, zb0, zb0+bzl))
H = nz*dx
ctrl0 = ("A 10 5\nB 50 0.002\nB 60 1\nW 10 %g\nF 50 2\nF 60 %g\nF 62 %g\nN 41 %g\nN 47 0.3\nM 10 NP\n"
         "P 10 0\nW 22 -9.81\n" % (Q, h0, h0, tend) + body)
RUNS = {"three": ([B_out, B_in, B_b], "C 11 21\nC 14 2\n", ctrl0 + "G 50 1\nG 54 1 1 %g\n" % Q),
        "one": ([one], "C 11 1\nC 14 2\n", ctrl0 + "G 50 1\nG 1 1\n")}

def control(c1, npr):
    xs = [B[0] for B in (B_in, B_b, B_out)]
    ys = [B[1] for B in (B_in, B_b, B_out)]
    x0, x1 = min(a.min() for a in xs), max(a.max() for a in xs)
    y0, y1 = min(a.min() for a in ys), max(a.max() for a in ys)
    return ("// canal with a 90 deg bend (t4_canal_bend.py)\n%sC 12 21\nC 13 21\nC 15 21\nC 16 3\n\nB 1 %g\n"
            "B 10 %g %g %g %g 0.0 %g\n\nM 10 %d\n\nR 20 1\n" % (c1, dx, x0, x1, y0, y1, H, npr))

which = os.environ.get("RUN", "both")
sys.stdout.reconfigure(line_buffering=True)
P = {}
for name, (bl, c1, ctrl) in RUNS.items():
    if which not in ("both", name):
        continue
    d = os.path.join(T6, "runs/t4_%s_np%d" % (name, np_))
    if not os.environ.get("NORUN"):
        rc = run4(d, control(c1, np_), ctrl.replace("NP", str(np_)), bl, np_=np_)
        print("%-5s rc=%d" % (name, rc))
    f = glob.glob(os.path.join(d, "REEF3D_*_6DOF", "REEF3D_6DOF_position_0.dat"))
    if f:
        rows = [l.split() for l in open(f[0]).readlines()[1:]]
        a = np.array([[float(v) for v in r] for r in rows if len(r) == 7])
        # up to 3 m before the outlet: the method cannot take the box out through the outflow boundary
        s_, _, _ = canal(a[:, 1], a[:, 2])
        out = np.nonzero(s_ > 0.5*np.pi*Rc + Lout - 3.0)[0]
        P[name] = a[:out[0]] if len(out) else a

for name, a in P.items():
    t, x, y, z, psi = a[:, 0], a[:, 1], a[:, 2], a[:, 3], a[:, 6]
    s, nn, ang = canal(x, y)
    sp = np.gradient(s, t)
    k_in, k_out = np.argmax(s > 0), np.argmax(s > 0.5*np.pi*Rc)
    print("%-5s t %.1f s: s %.2f m (%d samples)  n: max out %.3f m, max in %.3f m (wall at %.2f m from the axis, "
          "half the box %.2f)  speed (t > 20 s) mean %.3f max %.3f m/s (U0 %.2f)  heading psi %.1f .. %.1f deg (the axis turns by -90)  "
          "heave %.4f .. %.4f m  bend entered at t %.1f s, left at t %.1f s"
          % (name, t[-1], s[-1], len(t), nn.max(), -nn.min(), 0.5*Wc, 0.5*bxl, sp[t > 20].mean() if (t > 20).any() else np.nan,
             sp.max(), U0, psi.min(), psi.max(), z.min(), z.max(), t[k_in] if s.max() > 0 else np.nan,
             t[k_out] if s.max() > 0.5*np.pi*Rc else np.nan))
if "one" in P and "three" in P:
    a, b = P["one"], P["three"]
    tt = np.linspace(0, min(a[-1, 0], b[-1, 0]), 400)
    d = np.hypot(np.interp(tt, a[:, 0], a[:, 1]) - np.interp(tt, b[:, 0], b[:, 1]),
                 np.interp(tt, a[:, 0], a[:, 2]) - np.interp(tt, b[:, 0], b[:, 2]))
    print("three vs one: the position max %.4f m (at t %.1f s), at the end %.4f m; heading at the end %.2f deg"
          % (d.max(), tt[np.argmax(d)], d[-1], np.interp(tt[-1], b[:, 0], b[:, 6]) - np.interp(tt[-1], a[:, 0], a[:, 6])))

if os.environ.get("PNG") and P:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    C = ["#2B4C7E", "#E07B39", "#4F8FCF", "#3A9A5B", "#6B6B6B"]
    fig, ax = plt.subplots(1, 2, figsize=(12, 5.2), gridspec_kw={"width_ratios": [1.1, 1]})
    for B, lab in ((B_in, "block 1"), (B_b, "block 2"), (B_out, "block 0")):
        X, Y = B[0], B[1]
        for k in (0, -1):
            ax[0].plot(X[:, k], Y[:, k], "-", color=C[4], lw=1.2)
        ax[0].plot(X[0, :], Y[0, :], "--", color=C[4], lw=0.7)
        ax[0].plot(X[::8, :].T, Y[::8, :].T, "-", color="#cccccc", lw=0.4)
    for m, (name, a) in enumerate(P.items()):
        ax[0].plot(a[:, 1], a[:, 2], "-", color=C[m], lw=1.6, label="%s block%s" % (name, "s" if name == "three" else ""))
        for tq in range(0, int(a[-1, 0]) + 1, 10):
            k = np.argmin(np.abs(a[:, 0] - tq))
            ax[0].plot(a[k, 1], a[k, 2], "o", color=C[m], ms=3)
        s, nn, ang = canal(a[:, 1], a[:, 2])
        ax[1].plot(s, nn, "-", color=C[m], lw=1.4, label=name)
    ax[0].set_aspect("equal"); ax[0].set_xlabel("x (m)"); ax[0].set_ylabel("y (m)")
    ax[0].legend(frameon=False, fontsize=8, loc="lower right"); ax[0].set_title("path of the centre (dots every 10 s)", fontsize=9)
    for v in (0, 0.5*np.pi*Rc):
        ax[1].axvline(v, color=C[4], lw=0.8, ls="--")
    ax[1].axhline(0.5*Wc - 0.5*bxl, color=C[4], lw=0.8, ls=":")
    ax[1].axhline(-0.5*Wc + 0.5*bxl, color=C[4], lw=0.8, ls=":")
    ax[1].set_xlabel("s along the axis (m), bend between the dashed lines"); ax[1].set_ylabel("n from the axis (m, + outer bank)")
    ax[1].set_title("dotted: the box touches the bank", fontsize=9)
    for A in ax:
        A.grid(alpha=0.3)
        for sd in ("top", "right"):
            A.spines[sd].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.environ["PNG"], dpi=150)
    print("wrote", os.environ["PNG"])
