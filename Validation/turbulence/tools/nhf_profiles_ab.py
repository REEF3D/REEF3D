# Architect: Hans Bihs
"""nhf_profiles_ab.py <out.png> <label>=<run dir> ... [--x 220]
NHFLOW 2D channel: u(z) and nu_t(z) at one x station from the last VTU of each run, one curve per run.
Labels starting with 'before' are drawn dashed."""
import sys, os, glob, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import r3read

def profile(run, xq):
    d = r3read.read_vtk(sorted(glob.glob(run + '/REEF3D_NHFLOW_VTU/*-000001.vtu'))[-1])
    P = d['points']; m = P[:, 1] == 0.0
    xs = np.unique(np.round(P[m, 0], 6)); x = xs[np.argmin(abs(xs - xq))]
    sel = m & np.isclose(P[:, 0], x); o = np.argsort(P[sel, 2]); z = P[sel, 2][o]
    return z - z.min(), d['velocity'][sel][o, 0], d['eddyv'][sel][o]

if __name__ == '__main__':
    args = sys.argv[1:]; xq = 220.0
    if '--x' in args:
        i = args.index('--x'); xq = float(args[i + 1]); del args[i:i + 2]
    out, runs = args[0], [a.split('=', 1) for a in args[1:]]
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.5))
    for n, (lab, run) in enumerate(runs):
        z, u, nt = profile(run, xq)
        ls = '--' if lab.startswith('before') else '-'
        ax[0].plot(u, z, ls, color=f'C{n % 2 * 3}', label=lab); ax[1].plot(nt, z, ls, color=f'C{n % 2 * 3}', label=lab)
    ax[0].set_xlabel('u [m/s]'); ax[1].set_xlabel(r'$\nu_t$ [m$^2$/s]')
    for a in ax: a.set_ylabel('z [m]'); a.grid(alpha=.3); a.legend()
    fig.suptitle(f'NHFLOW 2D channel, x = {xq:g} m'); fig.tight_layout(); fig.savefig(out, dpi=110)
