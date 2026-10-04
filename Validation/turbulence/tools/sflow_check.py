# Architect: Hans Bihs
"""Compare SFLOW k, eps/omega, nu_t with the local depth-averaged equilibrium (Rastogi & Rodi):
   cf = g n^2/h^(1/3), n = ks^(1/6)/20 (as sflow_rough_manning), u* = sqrt(cf)|U|
   k = u*^2/(ceg sqrt(cmu) cf^(1/4)), eps = u*^3/(sqrt(cf) h), omega = eps/(cmu k), nu_t = cmu k^2/eps
"""
import sys, glob, numpy as np
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import r3read
def check(run, ks=0.01, ceg=2.7, cmu=0.09, g=9.81, model='ke', xr=(100, 900)):
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
    for run in sys.argv[1:]:
        model = 'kw' if run.endswith('_kw') else ('parab' if 'parab' in run else 'ke')
        r = check(run, model=model)
        if r is None: print(run, 'no output'); continue
        print(run)
        for k, v in r.items():
            if k.startswith('_'): continue
            print(f"   {k:22s} mean {np.nanmean(v):8.4f}  min {np.nanmin(v):8.4f}  max {np.nanmax(v):8.4f}")
