# P5b T6: the transverse bed slope of a river bend (SFLOW, curvilinear grid, S 86 and S 95), geometry and flow of the
# LFM flume of Delft Hydraulics (Struiksma et al. 1985): 180 deg, centre-line radius R 4.25 m, width B 1.7 m, depth
# h 0.2 m, discharge Q 0.17 m3/s, d50 0.78 mm, Chezy C 26.4 (Keulegan ks from it), bed slope 1.8 permille; straight
# reaches of 11 m upstream and 6 m downstream. The first 2 m of the inflow reach keep their bed and feed the bed
# load (S 71). After the bed has adapted: the transverse slope dz/dn of the cross-sections against the axisymmetric
# equilibrium of the model,
#     dz/dn = f(theta) alpha_I h / r,   f = A_sh theta^B_sh (d/h)^C_sh,   alpha_I = E_spir 2/kappa^2 (1 - sqrt(g)/(kappa C)),
# with the local depth, Shields parameter and Chezy coefficient of the run (the transverse bed-load flux vanishes).
# environment: XCTRL extra ctrl lines (;), ESPIR (S 96, 0.4), NOSPIRAL=1 (S 95 0), NOSLOPE=1 (S 86 0), MF (S 35, 100),
# DS (cell length, 0.1), NORUN=1
import sys, numpy as np
from p5lib import *
steps = int(sys.argv[1]) if len(sys.argv) > 1 else 16000
R, B, h0, Q, d50, Cz, slope = 4.25, 1.7, 0.2, 0.17, 0.00078, 26.4, 0.0018
L1, L2 = 11.0, 6.0
ds = float(os.environ.get("DS", "0.1"))
n1, n2, nb = int(round(L1/ds)), int(round(L2/ds)), int(round(np.pi*R/ds))
ny = int(round(B/ds)) | 1
nx = n1 + nb + n2
g, kappa, A, Bs, Cs = 9.81, 0.4, 9.0, 0.5, 0.3
espir = float(os.environ.get("ESPIR", "0.4"))     # E_spir 1 gives a point bar above the water (dry inner bank)
ks = 12.0*h0/np.exp(Cz/np.sqrt(g)/2.5)            # Keulegan: C = sqrt(g) 2.5 ln(12h/ks)
Delta = 2650.0/998.2 - 1.0

# nodes: s along the centre line, n to the left (towards the centre of the bend)
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
L = s_n[-1]

bc = "C 11 1\nC 12 21\nC 13 21\nC 14 2\nC 15 21\nC 16 3\n"
ctl = control(nx, ny, 1, nx*ds, ny*ds, 0.5, bc=bc)
ctrl = ("A 10 2\nA 210 3\nA 211 4\nA 220 0\nA 240 1\nB 50 %g\nB 60 1\nW 10 %g\nF 60 %g\nF 62 %g\nN 47 0.3\nM 10 1\nP 10 0\nW 22 -9.81\n"
        % (ks, Q, h0, h0)
        + "S 10 1\nS 11 2\nS 16 2\nS 20 %g\nS 21 %g\nS 15 3\nS 17 1.0\nS 35 %s\nS 41 1\nS 43 4000\nS 44 1\nS 31 1\nS 100 0\nS 71 2.0\n"
        % (d50, ks/d50, os.environ.get("MF", "100"))
        + "S 86 %d\nS 87 %g\nS 88 %g\nS 89 %g\nS 95 %d\nS 96 %g\n" % (0 if os.environ.get("NOSLOPE") else 1, A, Bs, Cs,
                                                                     0 if os.environ.get("NOSPIRAL") else 1, espir)
        + "G 50 2\n" + os.environ.get("XCTRL", "").replace(";", "\n") + "\n")
tag = os.environ.get("TAG", "")
d = os.path.join(S5, "runs/t6_bend%s" % tag)
rc = 0
if not os.environ.get("NORUN"):
    rc = run(d, ctl, ctrl, grid=(X, Y, Z), steps=steps)
st = final_state(d)
print("rc=%d t=%.1f s  ks %.4f m" % (rc, st["simtime"], ks),
      log_lines(d, "!!!")[:2])

# cell fields (i along s, j across)
bed = st["bed"].reshape(nx, ny)
WL = st["WL"].reshape(nx, ny)
z0 = 0.25*(Z[:-1, :-1] + Z[1:, :-1] + Z[:-1, 1:] + Z[1:, 1:])
P = st["P"].reshape(nx-1, ny)                                    # the inner faces i+1/2 (normal along s)
Qf = st["Q"].reshape(nx, ny-1)                                   # the inner faces j+1/2 (normal across)
Pe = np.vstack([P[:1], P, P[-1:]])
Qe = np.hstack([np.zeros((nx, 1)), Qf, np.zeros((nx, 1))])
U = 0.5*(Pe[:-1] + Pe[1:])
V = 0.5*(Qe[:, :-1] + Qe[:, 1:])
Umag = np.hypot(U, V)
nc = 0.5*(n_n[:-1] + n_n[1:])
sc = 0.5*(s_n[:-1] + s_n[1:])

def theory(i):
    """dz/dn of the axisymmetric equilibrium in section i (cell values), r the radius of the cells"""
    hh = np.maximum(WL[i], 1e-4)
    Cl = np.sqrt(g)*2.5*np.log(np.maximum(12*hh/ks, np.e))
    th = (np.sqrt(g)*Umag[i]/Cl)**2/(Delta*g*d50)
    f = A*np.sqrt(th)**(2*Bs)*(d50/hh)**Cs
    al = espir*2/kappa**2*(1 - np.sqrt(g)/(kappa*Cl))
    rr = R - nc
    return f*al*hh/rr, th, Cl

print("cross-sections (bend angle; bed change at the outer and inner bank; dz/dn of the run (central 3/5) vs the equilibrium):")
rows = []
core = slice(ny//5, ny - ny//5)
for ang in [0, 30, 45, 60, 90, 120, 150, 180]:
    i = n1 + min(nb-1, int(ang/180.0*nb))
    db = bed[i] - z0[i]
    ds_ = np.polyfit(nc[core], bed[i][core], 1)[0]
    th_, th, Cl = theory(i)
    rows.append((ang, ds_, th_[core].mean()))
    print("  %3d deg: outer %+.4f inner %+.4f m   dz/dn run %.4f  equilibrium %.4f  (ratio %.2f; theta %.2f C %.1f h %.3f)"
          % (ang, db[0], db[-1], ds_, th_[core].mean(), ds_/max(th_[core].mean(), 1e-12), th[core].mean(), Cl[core].mean(), WL[i][core].mean()))
dev = [a for a in rows if a[0] >= 90]
print("fully developed (90-180 deg): mean ratio run/equilibrium %.2f" % np.mean([a[1]/a[2] for a in dev]))
np.savez(os.path.join(d, "t6.npz"), bed=bed, WL=WL, U=U, V=V, z0=z0, sc=sc, nc=nc, n1=n1, nb=nb)
