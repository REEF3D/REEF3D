# Architect: Hans Bihs
"""fs_profile.py <run> [<run> ...]: CFD 2D channel (6 m), vertical profiles of nu_t and k at x = 5 m from the
regression dump; prints nu_t/(kappa u* h) against z/h near the free surface and band diagnostics."""
import sys, os, glob, re, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import r3read
H, U, KS = 0.24, 0.5, 0.001
US = U / (2.5 * np.log(11 * H / KS))
def grid(run):
    t = open(run + '/control.txt').read()
    dx = float(re.search(r'^B 1 ([\d.]+)', t, re.M).group(1))
    b = [float(v) for v in re.search(r'^B 10 (.*)$', t, re.M).group(1).split()]
    return dx, int(round((b[1] - b[0]) / dx)), int(round((b[5] - b[4]) / dx))
def profile(run, x=5.0):
    dx, nx, nz = grid(run)
    s = r3read.read_state(sorted(glob.glob(run + '/reg/state_*_r0.bin'))[-1])
    i = int(x / dx); sl = slice(i * nz, (i + 1) * nz)
    z = (np.arange(nz) + 0.5) * dx
    return dict(z=z, nut=s['eddyv'][sl], k=s['kin'][sl], phi=s['phi'][sl], dx=dx)
if __name__ == '__main__':
    for run in sys.argv[1:]:
        p = profile(run); zn = p['z'] / H; nn = p['nut'] / (0.4 * US * H)
        band = np.abs(p['phi']) < 1.6 * p['dx']
        water = p['phi'] > 0
        print(f"{run}: nut/(kappa u* h) at z/h:")
        for zz in (0.5, 0.7, 0.8, 0.9, 0.95, 1.0, 1.05):
            print(f"   z/h={zz:4.2f}  {np.interp(zz, zn, nn):7.4f}")
        print(f"   max in water {nn[water].max():.4f} at z/h={zn[water][np.argmax(nn[water])]:.3f};  max in band {nn[band].max():.4f};  max in air {nn[~water].max():.4f}")
