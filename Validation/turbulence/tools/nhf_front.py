# Architect: Hans Bihs
"""nhf_front.py <run> ... [--nx NX --nz NZ --x0 X0]: NHFLOW 2D (x-z) wetting front, from the regression-dump
states written every REEF3D_REGRESSION_EVERY steps (<run>/reg/state_*_r0.bin, cell values, single rank).
For every state: the shoreline (last column with h > 1 cm), the largest nu_t over all wet cells and over the
shallow columns (1 cm < h < 10 cm), and the largest turbulence length scale relative to the depth,
l/h = nu_t/(cmu^0.25 sqrt(k) h), in cells with k > 1e-8. Prints the maxima over time and when they occur."""
import sys, glob, struct, argparse, numpy as np

def read_state(fn):
    b = open(fn, 'rb').read()
    assert b[:8] == b'R3DREG01'
    pos = 8; rank, size, count = struct.unpack_from('<iii', b, pos); pos += 12
    t, = struct.unpack_from('<d', b, pos); pos += 8
    nf, = struct.unpack_from('<i', b, pos); pos += 4
    f = {}
    for _ in range(nf):
        name = b[pos:pos + 16].split(b'\0')[0].decode(); pos += 16
        n, = struct.unpack_from('<q', b, pos); pos += 8
        f[name] = np.frombuffer(b, '<f8', n, pos); pos += 8 * n
    return t, f

def front(run, nx, nz, x0, dx):
    cmu = 0.09
    out = []
    for fn in sorted(glob.glob(run + '/reg/state_*_r0.bin')):
        t, f = read_state(fn)
        if 'EV' not in f: continue
        ev = f['EV'].reshape(nx, nz); k = f['KIN'].reshape(nx, nz); h = f['WL'][:nx]
        wet = h > 0.01
        sh = wet & (h < 0.1)
        xs = x0 + (np.nonzero(wet)[0].max() + 0.5) * dx if wet.any() else np.nan
        ev_wet = ev[wet].max() if wet.any() else 0.0
        ev_sh = ev[sh].max() if sh.any() else 0.0
        l = np.where(k > 1e-8, ev / (cmu ** 0.25 * np.sqrt(np.maximum(k, 1e-30))), 0.0) / np.maximum(h, 1e-6)[:, None]
        l_sh = l[sh].max() if sh.any() else 0.0
        out.append((t, xs, ev_wet, ev_sh, l_sh))
    return np.array(out)

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('runs', nargs='+')
    ap.add_argument('--nx', type=int, default=400); ap.add_argument('--nz', type=int, default=10)
    ap.add_argument('--x0', type=float, default=0.0); ap.add_argument('--dx', type=float, default=0.1)
    a = ap.parse_args()
    for run in a.runs:
        r = front(run, a.nx, a.nz, a.x0, a.dx)
        print(run, f'({len(r)} states, t = {r[0,0]:.1f}-{r[-1,0]:.1f} s)')
        for c, name in ((1, 'shoreline x (m)'), (2, 'max nu_t wet'), (3, 'max nu_t, 1-10 cm'), (4, 'max l/h, 1-10 cm')):
            m = np.nanargmax(r[:, c])
            print(f'   {name:20s} max {r[m, c]:.4g} at t = {r[m, 0]:.1f} s   (median over time {np.nanmedian(r[:, c]):.4g})')
