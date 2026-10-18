# Architect: Hans Bihs
"""sflow_bank.py [--nx 200] [--dx 2.0] <run> ...: SFLOW 2D channel with a sloping bank (shoreline along x), from the
regression dump (last state; cells i-major, NX cells in x). Profiles across the width of the velocity u (P,
averaged to the cell centres), the depth h and nu_t over x0 < x < x1, and the ratio to the local Manning velocity
U_eq = h^(2/3) S^(1/2)/n (n = ks^(1/6)/20, S from the surface slope in x). Without lateral momentum exchange
U/U_eq = 1 everywhere; the diffusion lowers it next to a no-slip boundary."""
import sys, os, glob, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import r3read

def profile(run, nx=200, dx=2.0, x0=250.0, x1=350.0, ks=0.01, hmin=0.02):
    s = r3read.read_state(sorted(glob.glob(run + '/reg/state_*_r0.bin'))[-1])
    ny = s['WL'].size // nx
    h = s['WL'].reshape(nx, ny); eta = s['eta'].reshape(nx, ny); nut = s['eddyv'].reshape(nx, ny)
    P = s['P'].reshape(nx - 1, ny); u = np.vstack([P[:1], 0.5 * (P[1:] + P[:-1]), P[-1:]])
    x = (np.arange(nx) + 0.5) * dx
    sel = (x > x0) & (x < x1)
    wet = h[sel, 0] > hmin
    S = -np.polyfit(x[sel][wet], eta[sel, 0][wet], 1)[0]
    n = ks ** (1 / 6) / 20.0
    out = []
    for j in range(ny):
        hh = h[sel, j].mean(); uu = u[sel, j].mean()
        ueq = hh ** (2 / 3) * np.sqrt(max(S, 1e-12)) / n if hh > hmin else np.nan
        out.append((j + 0.5, hh, uu, uu / ueq if hh > hmin else np.nan, nut[sel, j].mean()))
    return S, np.array(out), s['_simtime']

if __name__ == '__main__':
    args = sys.argv[1:]; nx, dx = 200, 2.0
    for opt in ('--nx', '--dx'):
        if opt in args:
            i = args.index(opt); v = float(args[i + 1]); del args[i:i + 2]
            if opt == '--nx': nx = int(v)
            else: dx = v
    for run in args:
        S, a, t = profile(run, nx, dx)
        print(f"{run}   t = {t:.0f} s   S = {S:.3e}")
        print("      y       h        u      u/U_eq     nu_t")
        for r in a:
            if r[0] >= 20.0:
                print(f"   {r[0]:6.2f}  {r[1]:6.3f}  {r[2]:7.4f}  {r[3]:7.3f}  {r[4]:9.3e}")
