# Architect: Hans Bihs
"""cfd_fields.py [--box x0 x1] <run> ...: CFD 2D channel (NX x 18 cells, dx = 0.02 m) from the regression dump
(last state). Prints min / mean / max of k and nu_t in the water (phi > 0), min k in the interface band
(|phi| < 2 dx) and in the air, and with --box the means of u, k, nu_t in the water for x0 < x < x1."""
import sys, os, glob, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import r3read
NZ, DX = 18, 0.02

def load(run):
    global NX
    s = r3read.read_state(sorted(glob.glob(run + '/reg/state_*_r0.bin'))[-1])
    NX = s['phi'].size // NZ
    f = {k: v.reshape(NX, NZ) for k, v in s.items() if not k.startswith('_') and getattr(v, 'size', 0) == NX * NZ}
    if getattr(s.get('u'), 'size', 0) == (NX - 1) * NZ:   # u on the faces between the cells
        u = s['u'].reshape(NX - 1, NZ); f['u'] = np.vstack([u, u[-1:]])
    f['_t'] = s['_simtime']
    return f

if __name__ == '__main__':
    args = sys.argv[1:]; box = None
    if '--box' in args:
        i = args.index('--box'); box = (float(args[i + 1]), float(args[i + 2])); del args[i:i + 3]
    for run in args:
        f = load(run); phi = f['phi']; x = (np.arange(NX) + 0.5) * DX
        wet = phi > 0; band = np.abs(phi) < 2 * DX; air = phi < -2 * DX
        k, nu = f['kin'], f['eddyv']
        print(f"{run}  t={f['_t']:.2f}")
        print(f"   water: k min {k[wet].min():.3e} mean {k[wet].mean():.3e} max {k[wet].max():.3e} | "
              f"nu_t min {nu[wet].min():.3e} mean {nu[wet].mean():.3e} max {nu[wet].max():.3e}")
        print(f"   interface band: k min {k[band].min():.3e}   air: k min {k[air].min() if air.any() else float('nan'):.3e}")
        if box:
            m = wet & ((x > box[0]) & (x < box[1]))[:, None]
            o = wet & ((x > 1.0) & (x < box[0] - 0.5))[:, None]
            print(f"   box {box[0]}-{box[1]} m: u {f['u'][m].mean():.4f} k {k[m].mean():.3e} nu_t {nu[m].mean():.3e}   "
                  f"upstream 1-{box[0]-0.5} m: u {f['u'][o].mean():.4f} k {k[o].mean():.3e} nu_t {nu[o].mean():.3e}")
