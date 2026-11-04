#!/usr/bin/env python3
"""LIP 11D test 1B (Delta Flume 1993, erosive): SEASTATE surfbeat + SFLOW or NHFLOW with waves in the sediment transport.
usage: mk_lip.py <dir> <sflow|nhflow> [key=value ...]"""
import sys, os, numpy as np
d, model = sys.argv[1], sys.argv[2]; extra=sys.argv[3:]
os.makedirs(d, exist_ok=True)
src=os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','02_lip11d_1b','refdata','LIP11D_1B_initial_profile.csv')
a=np.loadtxt(src,delimiter=',',comments='#',skiprows=3)
x=a[:,0]; z=a[:,3]
L=float(x.max()); W=1.0; dx=1.0
with open(os.path.join(d,'geo.dat'),'w') as f:
    for xi,zi in zip(x,z):
        for y in (-1.0,0.5,2.0):
            f.write('%.3f %.3f %.4f\n'%(xi,y,zi))
if model=='sflow':
    ctl=['C 11 21','C 12 21','C 13 21','C 14 21','C 15 21','C 16 21','',
         'B 1 %.2f'%dx,'B 10 0.0 %.1f 0.0 %.1f -0.5 7.0'%(L,W),'','G 10 1','G 9 1','G 15 2','','M 10 1','M 20 2']
else:
    ctl=['C 11 21','C 12 21','C 13 21','C 14 21','C 15 21','C 16 21','',
         'B 2 %d 1 5'%int(L/dx),'B 10 0.0 %.1f 0.0 %.1f -0.5 7.0'%(L,W),'','G 10 1','G 9 1','G 15 2','','M 10 1','M 20 2']
open(os.path.join(d,'control.txt'),'w').write('\n'.join(ctl)+'\n')
wave=['A 700 1','A 701 30','A 702 0.08 0.6','A 703 24','A 705 0.01','A 710 0','A 711 1','A 712 1 0 2 2',
      'A 740 2','A 741 1.0 0.45','A 748 1','A 750 1','A 751 1','A 752 300.0','A 753 2',
      'A 770 1','A 772 3600.0','A 773 1','A 774 1','A 775 2',
      'B 85 2','B 88 3.3','B 93 1.4 5.0','B 130 0','B 131 0.0']
sed=['S 10 1','S 11 6','S 20 0.00022','S 21 3.0','S 22 2650.0','S 24 0.4','S 15 3','S 17 1.0','S 35 10.0',
     'S 41 2','S 45 600.0','S 44 1','S 50 4','S 110 1','S 116 1','S 117 0.1','S 118 1','S 119 0','S 90 1']
flow=['F 60 4.1','W 22 -9.81','M 10 1','P 10 1','P 30 1080.0','P 76 1','P 12 5000','P 182 1080.0']
tend=600.0+18*3600/10.0
if model=='sflow':
    ctrl=['A 10 2','A 218 1','B 50 0.01','S 16 2','S 120 2','N 41 %.1f'%tend,'N 47 0.4']+wave+sed+flow
else:
    ctrl=['A 10 5','A 534 1','A 535 0.55','A 536 2.0','A 560 1','B 50 0.01','S 28 2','S 16 1',
          'N 41 %.1f'%tend,'N 47 0.4','P 180 1','P 182 1080.0']+wave+sed+flow
keys={}
for l in ctrl: keys[' '.join(l.split()[:2])]=l
for e in extra:
    k,v=e.split('='); k=k.strip(); v=v.strip()
    if v=='del': keys.pop(k,None)
    else: keys[k]=k+' '+v
open(os.path.join(d,'ctrl.txt'),'w').write('\n'.join(keys.values())+'\n')
print(d,model,'L',L,'tend',tend)
