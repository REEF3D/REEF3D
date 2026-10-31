import sys,glob; import os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from vtu import read_vtu
from sed_eval import bedprof
import numpy as np
def nh_flux(run, xs):
    d=read_vtu(sorted(glob.glob(run+'/REEF3D_NHFLOW_VTU/*-000001.vtu'))[-1]); P=d['points']; m=P[:,1]<P[:,1].min()+1e-6
    out=[]
    for x in xs:
        mm=m&(abs(P[:,0]-x)<1e-3); o=np.argsort(P[mm,2]); z=P[mm,2][o]; u=d['velocity'][mm][o][:,0]; c=d['ST_conc'][mm][o]
        out.append((x, np.trapezoid(u*c,z), np.trapezoid(u,z), np.trapezoid(c,z)/(z[-1]-z[0]), z[-1]-z[0]))
    return out
if __name__=='__main__':
    run=sys.argv[1]; xs=[float(v) for v in sys.argv[2:]]
    x,zb,r,n=bedprof(run)
    if 'NHFLOW' in str(glob.glob(run+'/REEF3D_*')):
        for xx,qs,q,cm,h in nh_flux(run,xs):
            i=np.argmin(abs(x-xx)); print('x %6.1f  qs %.3e m2/s (%.4f kg/s/m)  q %.4f  Cbar %.3e  h %.3f  qbe %.3e (%.4f kg/s/m)  ub %.3f'%(xx,qs,qs*2650,q,cm,h,r['ST_qbe'][i],r['ST_qbe'][i]*2650,r['velocity'][i][0]))
    else:
        for xx in xs:
            i=np.argmin(abs(x-xx)); print('x %6.1f qbe %.3e (%.4f kg/s/m) cb %.3e cbe %.3e qbs %.3e (%.4f kg/s/m) u %.3f'%(xx,r['ST_qbe'][i],r['ST_qbe'][i]*2650,r['ST_cb'][i],r['ST_cbe'][i],r['ST_qbs'][i],r['ST_qbs'][i]*2650,r['velocity'][i][0]))
