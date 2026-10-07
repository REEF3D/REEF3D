# Architect: Hans Bihs
"""nhf_bank.py <run> ... [--x0 X0 --x1 X1 --y0 Y0]: NHFLOW 3D channel with a sloping bank, last VTU (vertex
values). Depth-averaged u, k and nu_t over the water column for each vertex row y >= y0 (default 20 m),
averaged over x0 < x < x1 (default 100-180 m), and the largest nu_t and nu_t/(k/omega or k^2/eps) ratio in
the wet region. Rows whose depth is below 1 cm are reported as dry."""
import sys, os, glob, argparse, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import r3read

def bank(run, x0, x1, y0):
    d = r3read.read_vtk(sorted(glob.glob(run + '/REEF3D_NHFLOW_VTU/*-000001.vtu'))[-1])
    P = d['points']; xs = np.round(P[:, 0], 4); ys = np.round(P[:, 1], 4)
    rows = []
    for y in np.unique(ys[ys >= y0]):
        acc = []
        for x in np.unique(xs[(xs > x0) & (xs < x1)]):
            c = (xs == x) & (ys == y)
            o = np.argsort(P[c, 2]); z = P[c, 2][o]; h = z[-1] - z[0]
            if h < 0.01: continue
            f = lambda a: np.trapezoid(a[c][o], z) / h
            acc.append((h, f(d['velocity'][:, 0]), f(d['kin']), f(d['eddyv'])))
        rows.append((y,) + (tuple(np.mean(acc, axis=0)) if acc else (0.0, np.nan, np.nan, np.nan)))
    return np.array(rows), d['eddyv'].max()

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('runs', nargs='+')
    ap.add_argument('--x0', type=float, default=100.0); ap.add_argument('--x1', type=float, default=180.0)
    ap.add_argument('--y0', type=float, default=20.0); a = ap.parse_args()
    for run in a.runs:
        r, numax = bank(run, a.x0, a.x1, a.y0)
        print(f"{run}   max nu_t = {numax:.3e}")
        for y, h, u, k, nu in r:
            print(f"   y = {y:5.1f}  h = {h:5.3f}" + ("   dry" if h == 0 else f"  u = {u:6.4f}  k = {k:.3e}  nu_t = {nu:.3e}"))
