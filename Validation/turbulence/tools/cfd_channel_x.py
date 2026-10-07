# Architect: Hans Bihs
"""cfd_channel_x.py <run> ...: CFD 2D open channel on any (also stretched) x grid, from the last VTU (vertex
values). Depth-averaged nu_t and k over the water (phi > 0) along x, divided by the equilibrium values
kappa u* h/6 and u*^2/(2 sqrt(cmu)) with u* = U/(2.5 ln(11 h/ks)); U = 0.5 m/s, h = 0.24 m, ks = 1 mm."""
import sys, os, glob, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import r3read

def profile(run, U=0.5, H=0.24, ks=0.001):
    d = r3read.read_vtk(sorted(glob.glob(run + '/REEF3D_CFD_VTU/*-000001.vtu'))[-1])
    P = d['points']; m = P[:, 1] == P[:, 1].min()
    us = U / (2.5 * np.log(11 * H / ks)); nu0 = 0.4 * us * H / 6; k0 = us ** 2 / 0.3 / 2
    out = []
    for x in np.unique(np.round(P[m, 0], 6)):
        c = m & np.isclose(P[:, 0], x) & (d['phi'] > 0)
        if c.sum() < 2: continue
        z = P[c, 2]; o = np.argsort(z); z = z[o]
        nu = d['eddyv'][c][o]; k = d['kin'][c][o]
        h = z.max() - z.min()
        out.append((x, np.trapezoid(nu, z) / h / nu0, np.trapezoid(k, z) / h / k0))
    return np.array(out)

if __name__ == '__main__':
    stations = (0.5, 2.0, 4.0, 6.0, 8.0, 10.0, 11.8)
    for run in sys.argv[1:]:
        a = profile(run)
        print(run)
        for xs in stations:
            r = a[np.argmin(abs(a[:, 0] - xs))]
            print(f"   x = {r[0]:6.2f}  nu_t/nu_eq = {r[1]:6.3f}  k/k_eq = {r[2]:6.3f}")
