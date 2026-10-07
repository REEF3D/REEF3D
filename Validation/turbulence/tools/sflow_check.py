# Architect: Hans Bihs
"""Compare SFLOW k, eps/omega, nu_t with the local depth-averaged equilibrium (Rastogi & Rodi):
   cf = g n^2/h^(1/3), n = ks^(1/6)/20 (as sflow_rough_manning), u* = sqrt(cf)|U|
   k = u*^2/(ceg sqrt(cmu) cf^(1/4)), eps = u*^3/(sqrt(cf) h), omega = eps/(cmu k), nu_t = cmu k^2/eps
ceg = A 264 from <run>/ctrl.txt, else the default (3.6 since the 2026-10 turbulence patch, 2.7 before;
--ceg <value> sets the default for runs made with an older build).
usage: sflow_check.py [--ceg 2.7] <run> ...
"""
import sys, glob, numpy as np
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import r3read
def a264(run, default):
    try:
        for line in open(os.path.join(run, 'ctrl.txt')):
            t = line.split()
            if len(t) >= 3 and t[0] == 'A' and t[1] == '264': return float(t[2])
    except OSError:
        pass
    return default

def check(run, ks=0.01, ceg=3.6, cmu=0.09, g=9.81, model='ke', xr=(100, 900)):
    fs = sorted(glob.glob(run + '/REEF3D_SFLOW_VTP_FSF/*-000001.vtp'))
    if not fs: return None
    d = r3read.read_vtk(fs[-1])
    x = d['points'][:, 0]; y = d['points'][:, 1]
    sel = (x > xr[0]) & (x < xr[1])
    U = np.linalg.norm(d['velocity'][:, :2], axis=1)[sel]
    h = d['waterlevel'][sel]
    n = ks ** (1 / 6) / 20.0; cf = g * n * n / h ** (1 / 3); us = np.sqrt(cf) * U
    keq = us ** 2 / (ceg * np.sqrt(cmu) * cf ** 0.25); eeq = us ** 3 / (np.sqrt(cf) * h)
    weq = eeq / (cmu * keq); nueq = cmu * keq ** 2 / eeq
    r = {'nut/nut_eq': d['eddyv'][sel] / nueq}
    if 'kin' in d: r['k/k_eq'] = d['kin'][sel] / keq
    if model == 'ke': r['eps/eps_eq'] = d['epsilon'][sel] / eeq
    if model == 'kw': r['omega/omega_eq'] = d.get('omega', d.get('epsilon'))[sel] / weq
    if model == 'parab': r['nut/(kappa/6 u* h)'] = d['eddyv'][sel] / (0.4 / 6 * us * h)
    r['_y'] = y[sel]
    return r
if __name__ == '__main__':
    args = sys.argv[1:]; ceg0 = 3.6
    if '--ceg' in args:
        i = args.index('--ceg'); ceg0 = float(args[i + 1]); del args[i:i + 2]
    for run in args:
        model = 'kw' if run.rstrip('/').endswith('_kw') else ('parab' if 'parab' in run else 'ke')
        r = check(run, ceg=a264(run, ceg0), model=model)
        if r is None: print(run, 'no output'); continue
        print(run)
        for k, v in r.items():
            if k.startswith('_'): continue
            print(f"   {k:22s} mean {np.nanmean(v):8.4f}  min {np.nanmin(v):8.4f}  max {np.nanmax(v):8.4f}")
