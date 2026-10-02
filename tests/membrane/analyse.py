#!/usr/bin/env python3
# Static closed bag (X 330): compares the membrane time series with the expected Darcy leakage and
# hydrostatic floor load.  usage: analyse.py REEF3D_NHFLOW_Membrane/REEF3D_NHFLOW_Membrane_0.dat R_n A_exposed
#   R_n        hydraulic resistance from membrane.dat [m/s]
#   A_exposed  membrane area carrying the pressure difference (walls below the free surface + floor) [m^2]
import sys
import numpy as np

d = np.loadtxt(sys.argv[1])
Rn = float(sys.argv[2])
Aexp = float(sys.argv[3])
t, ein, eout, dh, Q, Fx, Fy, Fz, Fzf, Fzh, urel, umax, vol = d.T[:13]
g = 9.81

print(f"{'t':>7} {'dh':>9} {'Q_leak':>10} {'Q_darcy':>10} {'Fz_floor':>10} {'-rho g dh A':>11} {'ratio':>6} {'max|U|':>9}")
for n in np.unique(np.searchsorted(t, np.linspace(t[0], t[-1], 11)).clip(0, len(t)-1)):
    print(f"{t[n]:7.2f} {dh[n]:9.5f} {Q[n]:10.3e} {g*dh[n]/Rn*Aexp:10.3e} {Fzf[n]:10.2f} {Fzh[n]:11.2f} {Fzf[n]/Fzh[n]:6.3f} {umax[n]:9.2e}")
print(f"water volume drift: {(vol[-1]-vol[0])/vol[0]:.2e}")
