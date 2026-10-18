# Architect: Hans Bihs
import sys,glob; import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import r3read, numpy as np
def profiles(run, xs, last=True):
    fs=sorted(glob.glob(run+'/REEF3D_NHFLOW_VTU/*-000001.vtu')); d=r3read.read_vtk(fs[-1])
    P=d['points']; m=P[:,1]==0.0
    out={}
    for x in xs:
        sel=m & (np.abs(P[:,0]-x)<1e-6); o=np.argsort(P[sel,2])
        out[x]={k:(d[k][sel][o] if d[k].ndim==1 else d[k][sel][o,0]) for k in ('velocity','eddyv','kin','epsilon','omega') if k in d}
        out[x]['z']=P[sel,2][o]
    return out
if __name__=='__main__' and False:
    run=sys.argv[1]; xs=[float(a) for a in sys.argv[2:]] or [3,30,150,285]
    pr=profiles(run,xs)
    for x,p in pr.items():
        print(f"x={x}: z  ",np.round(p['z'][::2],2))
        for k in ('velocity','kin','epsilon','omega','eddyv'):
            if k in p: print("   ",k,np.round(p[k][::2],5))

def summary(run, x=45.0, H=2.0, ks=0.01):
    p=profiles(run,[x])[x]; z=p['z']; u=p['velocity']; nu=p['eddyv']
    m=(z>0.05*H)&(z<0.3*H)
    A=np.vstack([np.log(30*z[m]/ks)/0.4,np.ones(m.sum())]).T
    us=np.linalg.lstsq(A,u[m],rcond=None)[0][0]
    zz=z[(z>0)&(z<H)]; nn=nu[(z>0)&(z<H)]
    nu_avg=np.trapezoid(nu,z)/H
    return dict(ustar_fit=us, nut_avg=nu_avg, nut_avg_theory=0.4*us*H/6, ratio=nu_avg/(0.4*us*H/6), U=np.trapezoid(u,z)/H, k_bed=p['kin'][1], k_bed_theory=us**2/0.3)

if __name__=='__main__':
    # nhf_profiles.py <run dir> [x ...]: depth-averaged nu_t against kappa u* h/6 (u* fitted to the log law)
    run=sys.argv[1]; xs=[float(a) for a in sys.argv[2:]] or [3,30,150,285]
    for x in xs:
        s=summary(run,x); print(f"x={x:6.1f} " + ' '.join(f'{k}={v:.4g}' for k,v in s.items()))
