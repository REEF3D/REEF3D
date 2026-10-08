# FW-H (permeable box) vs Curle compact dipole from the body force, and vs pressure probes
#   Curle: p(x,t) = -1/(4 pi) [ rhat.F/r^2 + rhat.dF/dt/(c r) ]_(t - r/c),  F = force of the fluid on the body
import math, sys, bisect
c0 = 1500.0
T0 = float(sys.argv[1]) if len(sys.argv)>1 else 25.0
T1 = float(sys.argv[2]) if len(sys.argv)>2 else 1e9
obs = [(0,50,0),(0,0,50),(50,0,0),(0,-50,0),(0,3,0),(0,0,3)]

def rows(fn, ncol):
    d=[]
    for l in open(fn):
        s=l.split()
        if len(s)<ncol: continue
        try: d.append([float(x) for x in s[:ncol]])
        except: pass
    return d

F = rows("REEF3D_CFD_6DOF/REEF3D_6DOF_forces_0.dat", 4)
# drop repeated times, keep the last value of each time
Fd={}
for r in F: Fd[r[0]]=r[1:4]
tF=sorted(Fd); FF=[Fd[t] for t in tF]
def Fat(t):
    k=bisect.bisect_left(tF,t)
    k=min(max(k,1),len(tF)-1)
    t0,t1=tF[k-1],tF[k]; w=(t-t0)/(t1-t0)
    v=[FF[k-1][i]*(1-w)+FF[k][i]*w for i in range(3)]
    dv=[(FF[k][i]-FF[k-1][i])/(t1-t0) for i in range(3)]
    return v,dv

def stats(a,b):
    ma=sum(a)/len(a); mb=sum(b)/len(b)
    a=[x-ma for x in a]; b=[x-mb for x in b]
    ra=math.sqrt(sum(x*x for x in a)/len(a)); rb=math.sqrt(sum(x*x for x in b)/len(b))
    cor=sum(x*y for x,y in zip(a,b))/len(a)/(ra*rb)
    rd=math.sqrt(sum((x-y)**2 for x,y in zip(a,b))/len(a))
    return ra,rb,cor,rd/rb

print(f"window {T0}..{min(T1,tF[-1]):.1f} s, fluctuations (means removed)")
print(f"{'observer':>14} {'rms FW-H':>10} {'rms Curle':>10} {'ratio':>6} {'corr':>6} {'rms diff/Curle':>14}")
for n,x in enumerate(obs):
    r=math.sqrt(sum(c*c for c in x)); rh=[c/r for c in x]
    fw=rows(f"REEF3D_CFD_Acoustics/REEF3D-CFD-FWH-Observer-{n+1}.dat",2)
    A=[];B=[]
    for t,p in fw:
        if t<T0 or t>T1 or t-r/c0>tF[-1]: continue
        v,dv=Fat(t-r/c0)
        cu=-(sum(rh[i]*v[i] for i in range(3))/(r*r) + sum(rh[i]*dv[i] for i in range(3))/(c0*r))/(4*math.pi)
        A.append(p); B.append(cu)
    if len(A)<10: print(x,"too few samples"); continue
    ra,rb,cor,rel=stats(A,B)
    print(f"{str(x):>14} {ra:10.3e} {rb:10.3e} {ra/rb:6.3f} {cor:6.3f} {rel:14.3f}")

print("near observers vs pressure probes P 64 (means removed)")
for n,(x,po) in enumerate([((0,3,0),5),((0,0,3),6)]):
    pr=rows(f"REEF3D_CFD_PressureProbe/REEF3D-CFD-Probe-Pressure-{n+1}.dat",2)
    fw=rows(f"REEF3D_CFD_Acoustics/REEF3D-CFD-FWH-Observer-{po}.dat",2)
    tp=[r[0] for r in pr]; pp=[r[1] for r in pr]
    A=[];B=[]
    for t,p in fw:
        if t<T0 or t>T1 or t>tp[-1]: continue
        k=bisect.bisect_left(tp,t); k=min(max(k,1),len(tp)-1)
        w=(t-tp[k-1])/(tp[k]-tp[k-1]); A.append(p); B.append(pp[k-1]*(1-w)+pp[k]*w)
    ra,rb,cor,rel=stats(A,B)
    print(f"{str(x):>14} rms FW-H {ra:.3e}  rms probe {rb:.3e}  ratio {ra/rb:.3f}  corr {cor:.3f}")
