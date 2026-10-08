# box momentum balance vs REEF3D body force
#   F_fwh = -(Lp + Lm + dM/dt)       dipole radiated by the FW-H integral (viscous stresses neglected)
#   F_box = -(Lp + Lm - V + dM/dt)   full momentum balance of the box
import math, sys, bisect
T0 = float(sys.argv[1]) if len(sys.argv)>1 else 25.0
T1 = float(sys.argv[2]) if len(sys.argv)>2 else 1e9
def rows(fn,n):
    d={}
    for l in open(fn):
        s=l.split()
        if len(s)<n: continue
        try: v=[float(x) for x in s[:n]]
        except: continue
        d[v[0]]=v[1:]
    t=sorted(d); return t,[d[x] for x in t]
tb,B=rows("REEF3D_CFD_Acoustics/REEF3D-CFD-FWH-Box.dat",13)
tf,F=rows("REEF3D_CFD_6DOF/REEF3D_6DOF_forces_0.dat",13)
def at(tt,V,t,i):
    k=bisect.bisect_left(tt,t); k=min(max(k,1),len(tt)-1)
    w=(t-tt[k-1])/(tt[k]-tt[k-1]); return V[k-1][i]*(1-w)+V[k][i]*w
cols={k:[] for k in ("Lp","Lm","V","dM","fwh","box","F","Fp")}
for k in range(1,len(tb)-1):
    t=tb[k]
    if t<T0 or t>T1 or t>tf[-1]: continue
    dM=[(B[k+1][9+q]-B[k-1][9+q])/(tb[k+1]-tb[k-1]) for q in range(3)]
    Lp=B[k][0:3]; Lm=B[k][3:6]; V=B[k][6:9]
    cols["Lp"].append([-x for x in Lp]); cols["Lm"].append([-x for x in Lm]); cols["V"].append(V); cols["dM"].append([-x for x in dM])
    cols["fwh"].append([-(Lp[q]+Lm[q]+dM[q]) for q in range(3)])
    cols["box"].append([-(Lp[q]+Lm[q]-V[q]+dM[q]) for q in range(3)])
    cols["F"].append([at(tf,F,t,q) for q in range(3)]); cols["Fp"].append([at(tf,F,t,6+q) for q in range(3)])
def ms(a):
    m=sum(a)/len(a); return m, math.sqrt(sum((x-m)**2 for x in a)/len(a))
def cor(a,b):
    ma=sum(a)/len(a); mb=sum(b)/len(b)
    return sum((x-ma)*(y-mb) for x,y in zip(a,b))/math.sqrt(sum((x-ma)**2 for x in a)*sum((y-mb)**2 for y in b))
n=len(cols["F"]); print(f"window {T0}..{min(T1,tb[-1]):.1f} s, {n} samples")
print("mean [N]   -Lp      -Lm      +V     -dM/dt |  F_fwh    F_box  | F REEF3D  (pressure part)")
for q,nm in enumerate("xyz"):
    m={k:ms([v[q] for v in cols[k]])[0] for k in cols}
    print(f"  F{nm} {m['Lp']:8.2f} {m['Lm']:8.2f} {m['V']:7.2f} {m['dM']:8.2f} | {m['fwh']:8.2f} {m['box']:8.2f} | {m['F']:8.2f} ({m['Fp']:.2f})")
print("fluct. rms  F_fwh   F_box   F      fwh/F  box/F  corr(fwh,F)  rms(V)")
for q,nm in enumerate("xyz"):
    a=[v[q] for v in cols["fwh"]]; b=[v[q] for v in cols["box"]]; c=[v[q] for v in cols["F"]]; vv=[v[q] for v in cols["V"]]
    ra,rb,rc,rv=ms(a)[1],ms(b)[1],ms(c)[1],ms(vv)[1]
    print(f"  F{nm} {ra:7.3f} {rb:7.3f} {rc:7.3f}  {ra/rc:6.3f} {rb/rc:6.3f}  {cor(a,c):6.3f}      {rv:6.3f}")
