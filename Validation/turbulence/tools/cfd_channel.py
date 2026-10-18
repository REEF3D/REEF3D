# Architect: Hans Bihs
"""CFD 2D open channel: k, eps/omega, nu_t from the regression dump (cell centres).
Compares with the equilibrium log-law profile u* = U/(2.5 ln(11H/ks)), nu_t = kappa u* z (1-z/H), k = u*^2/sqrt(cmu)(1-z/H)."""
import sys, glob, numpy as np
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import r3read
NX, NZ, DX = 600, 18, 0.02
def load(run):
    fs = sorted(glob.glob(run + '/reg/state_*_r0.bin'))
    s = r3read.read_state(fs[-1])
    f = {k: v.reshape(NX, NZ) for k, v in s.items() if not k.startswith('_') and getattr(v, 'size', 0) == NX * NZ}
    f['_t'] = s['_simtime']; return f
def column(f, xi, H=0.24):
    wet = f['phi'][xi] > 0
    z = (np.arange(NZ) + 0.5) * DX
    nu = f['eddyv'][xi][wet]; k = f['kin'][xi][wet]
    return dict(nut_avg=np.sum(nu * DX) / H, k_avg=np.sum(k * DX) / H, k_min_wet=k.min(), eps_min=f['eps'][xi][wet].min())
def theory(U=0.5, H=0.24, ks=0.001):
    us = U / (2.5 * np.log(11 * H / ks)); return dict(ustar=us, nut_avg=0.4 * us * H / 6, k_avg=us ** 2 / 0.3 / 2)
if __name__ == '__main__':
    th = theory(); print('theory', {k: round(v, 6) for k, v in th.items()})
    for run in sys.argv[1:]:
        f = load(run); print(run, 't=', round(f['_t'], 2))
        for x in (0.1, 0.5, 1, 2, 4, 6, 8, 11.8):
            c = column(f, int(x / DX)); print(f"   x={x:5.1f}  nut_avg/th={c['nut_avg']/th['nut_avg']:6.3f}  k_avg/th={c['k_avg']/th['k_avg']:6.3f}  min k={c['k_min_wet']:.2e}  min eps/omega={c['eps_min']:.3g}")
