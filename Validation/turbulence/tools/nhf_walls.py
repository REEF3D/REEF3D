# Architect: Hans Bihs
"""nhf_walls.py <run> ...: NHFLOW 3D channel with side walls. Profiles across the width (y) of U, k and nu_t,
averaged over x0 < x < x1 and 0.3 h < z < 0.7 h, from the last VTU. Prints the profile and the ratio of the
wall-adjacent column to the centre column."""
import sys, os, glob, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import r3read

def profiles(run, x0=80.0, x1=110.0):
    d = r3read.read_vtk(sorted(glob.glob(run + '/REEF3D_NHFLOW_VTU/*-000001.vtu'))[-1])
    P = d['points']; x, y, z = P[:, 0], P[:, 1], P[:, 2]
    zb = z.min(); h = np.full_like(z, np.nan)
    # water depth per (x, y) column from the highest vertex
    key = np.round(x, 6) * 1.0e6 + np.round(y, 6)
    for kk in np.unique(key):
        m = key == kk; h[m] = z[m].max() - zb
    rel = (z - zb) / h
    sel = (x > x0) & (x < x1) & (rel > 0.3) & (rel < 0.7)
    ys = np.unique(np.round(y[sel], 6))
    out = []
    for yy in ys:
        m = sel & np.isclose(y, yy)
        out.append((yy, d['velocity'][m, 0].mean(), d['kin'][m].mean(), d['eddyv'][m].mean()))
    return np.array(out)

if __name__ == '__main__':
    for run in sys.argv[1:]:
        a = profiles(run)
        c = a[len(a) // 2]
        print(run)
        print("      y        U         k          nu_t")
        for r in a:
            print(f"   {r[0]:6.3f}  {r[1]:7.4f}  {r[2]:9.3e}  {r[3]:9.3e}")
        print(f"   wall/centre:  U {a[0,1]/c[1]:.3f}  k {a[0,2]/c[2]:.3f}  nu_t {a[0,3]/c[3]:.3f}"
              f"   (other wall: U {a[-1,1]/c[1]:.3f}  k {a[-1,2]/c[2]:.3f}  nu_t {a[-1,3]/c[3]:.3f})")
