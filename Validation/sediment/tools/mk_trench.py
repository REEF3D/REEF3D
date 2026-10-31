#!/usr/bin/env python3
"""van Rijn (1980/1986) trench: write control.txt, ctrl.txt and geo.dat for SFLOW or NHFLOW.
usage: mk_trench.py <dir> <sflow|nhflow> <test 1|2|3> [key=value ...]  (extra ctrl lines 'S 35=10')"""
import sys, os
d, model, test = sys.argv[1], sys.argv[2], int(sys.argv[3])
extra = [e for e in sys.argv[4:] if not e.startswith('S 35=')]
S35 = float(([e.split('=')[1] for e in sys.argv[4:] if e.startswith('S 35=')] or ['10'])[0])
os.makedirs(d, exist_ok=True)
# geometry (station-1 coordinates, z below the flat bed): van Rijn 1986b Fig. 16 (ANUGA digitisation)
geo = {0: [(0.0,0.0),(6.0,0.0)],
       1: [(0.0,0.0),(1.0,0.0),(1.5,0.175),(5.5,0.175),(6.0,0.0)],
       2: [(0.0,0.0),(1.2,0.170),(3.8,0.170),(5.0,0.0)],
       3: [(0.0,0.0),(1.5,0.150),(4.5,0.150),(6.0,0.0)]}[test]
x0, x1 = -20.0, 20.0           # 20 m approach, 14 m beyond the trench
zb = 0.30                      # flat bed level (sand layer below)
h0, U0 = 0.39, 0.51
W = 0.5
q = h0*U0
dx = 0.1
nx = int(round((x1-x0)/dx))
def zbed(x):
    if x<=geo[0][0] or x>=geo[-1][0]: return zb
    for (xa,za),(xb,zb_) in zip(geo[:-1],geo[1:]):
        if xa<=x<=xb:
            return zb - (za + (zb_-za)*(x-xa)/(xb-xa) if xb>xa else za)
    return zb
# geo.dat: profile points every dx/2 at three y positions
with open(os.path.join(d,'geo.dat'),'w') as f:
    n=int(round((x1-x0)/(dx/2)))
    for i in range(n+1):
        x=x0+i*dx/2
        for y in (-0.5, 0.25, 1.0):
            f.write('%.4f %.4f %.5f\n'%(x,y,zbed(x)))
if model=='sflow':
    ctl = ['C 11 1','C 12 3','C 13 3','C 14 2','C 15 21','C 16 3','',
           'B 1 %.3f'%dx,'B 10 %.2f %.2f 0.0 %.2f 0.0 1.0'%(x0,x1,W),'',
           'G 10 1','G 9 1','G 15 2','','M 10 1','M 20 2']
else:
    ctl = ['C 11 1','C 12 3','C 13 3','C 14 2','C 15 21','C 16 3','',
           'B 2 %d 1 20'%nx,'B 10 %.2f %.2f 0.0 %.2f 0.0 %.2f'%(x0,x1,W,zb+h0),'',
           'G 10 1','G 9 1','G 15 2','','M 10 1','M 20 2']
open(os.path.join(d,'control.txt'),'w').write('\n'.join(ctl)+'\n')
rho_s=2650.0
qb_in = 0.01/rho_s                # measured bed load supply, m^2/s
qs_in = 0.03/rho_s                # measured suspended supply, m^2/s
common = ['S 10 1','S 11 1','S 15 3','S 17 1.0','S 20 0.00016','S 21 3.75','S 22 2650.0','S 24 0.4',
          'S 25 0','S 23 0.013','S 35 %.1f'%S35,'S 41 2','S 45 300.0','S 44 1','S 50 2',
          'S 129 1','S 130 %.4e'%qb_in,'S 90 1',
          'B 60 1','W 10 %.5f'%(q*W),'W 22 -9.81','F 60 %.3f'%(zb+h0),'F 62 %.3f'%(zb+h0),
          'N 41 %.1f'%(300.0+15*3600/max(S35,1e-9)),'M 10 1','P 10 1','P 30 %.1f'%(5400.0/max(S35,1e-9)),'P 76 1','P 26 1','P 12 2000','P 180 1','P 182 %.1f'%(5400.0/max(S35,1e-9))]
if model=='sflow':
    ctrl = ['A 10 2','A 220 0','A 226 1','A 228 0.1','A 218 1','B 50 0.025','S 16 2','S 120 1','S 124 %.4e'%(qs_in/q),'N 47 0.3','F 50 2','I 11 1']+common
else:
    ctrl = ['A 10 5','A 520 0','A 534 1','A 535 1.0','A 536 0.1','A 560 1','A 567 1','B 50 0.025','S 28 2','S 16 1','S 12 1','S 62 1',
            'S 124 %.4e'%(qs_in/q),'N 47 0.5']+common
keys={}
for l in ctrl: keys[' '.join(l.split()[:2])]=l
for e in extra:
    k,v=e.split('='); keys[k.strip()]=k.strip()+' '+v.strip()
    if v.strip()=='del': del keys[k.strip()]
open(os.path.join(d,'ctrl.txt'),'w').write('\n'.join(keys.values())+'\n')
print(d, model, 'test',test,'nx',nx,'q',q,'qb_in',qb_in,'C_in',qs_in/q)
