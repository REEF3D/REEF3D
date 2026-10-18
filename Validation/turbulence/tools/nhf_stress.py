# Architect: Hans Bihs
"""nhf_stress.py <run> ...: NHFLOW 2D channel, momentum balance of the developed flow.
Turbulent shear stress (nu+nu_t) du/dz against the driving stress g S (h - z), S = -d(eta)/dx
from the free surface between x0 and x1; ratio averaged over 0.15 h < z < 0.6 h."""
import sys, os, glob, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import r3read
def balance(run, x0=150.0, x1=280.0, nu=1.0e-6, g=9.81):
    d = r3read.read_vtk(sorted(glob.glob(run + '/REEF3D_NHFLOW_VTU/*-000001.vtu'))[-1])
    P = d['points']; m = P[:, 1] == 0.0
    xs = np.unique(np.round(P[m, 0], 6)); top = []
    for x in xs:
        sel = m & np.isclose(P[:, 0], x); top.append(P[sel, 2].max())
    top = np.array(top); fit = (xs >= x0) & (xs <= x1)
    S = -np.polyfit(xs[fit], top[fit], 1)[0]
    ratios = []
    for x in xs[fit][::4]:
        sel = m & np.isclose(P[:, 0], x); o = np.argsort(P[sel, 2])
        z = P[sel, 2][o]; u = d['velocity'][sel][o, 0]; nt = d['eddyv'][sel][o]
        h = z.max() - z.min(); zc = 0.5 * (z[1:] + z[:-1]) - z.min()
        tau = (nu + 0.5 * (nt[1:] + nt[:-1])) * np.diff(u) / np.diff(z)
        drv = g * S * (h - zc); w = (zc > 0.15 * h) & (zc < 0.6 * h)
        ratios.append(np.mean(tau[w] / drv[w]))
    return S, np.mean(ratios), np.std(ratios)
if __name__ == '__main__':
    for r in sys.argv[1:]:
        S, rm, rs = balance(r); print(f"{r:50s} S={S:.3e}  (nu+nu_t)du/dz / gS(h-z) = {rm:.3f} +- {rs:.3f}")
