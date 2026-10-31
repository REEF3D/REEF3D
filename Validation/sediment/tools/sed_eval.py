import sys,glob,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from vtp import read_vtp
REF=os.path.join(os.path.dirname(os.path.abspath(__file__)),'..')
import numpy as np
def bedprof(run, which=-1):
    sub='NHFLOW' if glob.glob(run+'/REEF3D_NHFLOW_VTP_BED') else 'SFLOW'
    fs=sorted(glob.glob(run+'/REEF3D_%s_VTP_BED/*-000001.vtp'%sub))
    if not fs: fs=sorted(glob.glob(run+'/REEF3D_%s_VTP_BED/*.vtp'%sub))
    b=read_vtp(fs[which]); P=b['?'].reshape(-1,3); y0=P[:,1].min(); m=abs(P[:,1]-y0)<1e-6; o=np.argsort(P[m,0])
    n=len(P); r={}
    for k in b:
        if k in ('?','connectivity','offsets'): continue
        a=b[k]; a=a.reshape(n,-1) if a.size%n==0 and a.size>n else a
        r[k]=a[m][o]
    return P[m,0][o], P[m,2][o], r, len(fs)
def trench_meas(test):
    a=np.loadtxt(os.path.join(REF,'01_vanrijn_trench','refdata','van_rijn_1986_fig16_measured.csv'),delimiter=',',comments='#')
    a=a[a[:,0]==test]; return a[:,1], a[:,2]
def trench_init(test):
    a=np.loadtxt(os.path.join(REF,'01_vanrijn_trench','refdata','van_rijn_1986_fig16_initial.csv'),delimiter=',',comments='#')
    a=a[a[:,0]==test]; return a[:,1], a[:,2]
def trench_score(run,test,zb=0.30):
    x,z,_,n=bedprof(run); depth=zb-z           # below the flat bed
    xm,zm=trench_meas(test); zi=np.interp(xm,x,depth)
    err=zi-zm; i=np.argmax(depth*((x>2)&(x<11)))
    return dict(rms=float(np.sqrt(np.mean(err**2))), emax=float(np.abs(err).max()), xdeep=float(x[i]), zdeep=float(depth[i]), nout=n)
def lip_meas():
    a=np.loadtxt(os.path.join(REF,'02_lip11d_1b','refdata','LIP11D_1C_initial_profile_approx_1B_end.csv'),delimiter=',',comments='#',skiprows=3)
    return a[:,0],a[:,3]
def lip_init():
    a=np.loadtxt(os.path.join(REF,'02_lip11d_1b','refdata','LIP11D_1B_initial_profile.csv'),delimiter=',',comments='#',skiprows=3)
    return a[:,0],a[:,3]
def lip_score(run,xr=(100,190)):
    x,z,_,n=bedprof(run); xm,zm=lip_meas(); xi,zi=lip_init()
    m=(xm>=xr[0])&(xm<=xr[1]); zc=np.interp(xm[m],x,z); z0=np.interp(xm[m],xi,zi)
    bss=1-np.mean((zc-zm[m])**2)/np.mean((z0-zm[m])**2)
    return dict(rms=float(np.sqrt(np.mean((zc-zm[m])**2))), bss=float(bss), nout=n)
if __name__=='__main__':
    kind,run=sys.argv[1],sys.argv[2]
    print(run, trench_score(run,int(sys.argv[3])) if kind=='trench' else lip_score(run))
