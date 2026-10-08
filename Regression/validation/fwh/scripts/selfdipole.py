# FW-H at the far observers vs the compact dipole of the run's own box momentum balance
#   F_fwh = -(Lp + Lm + dM/dt) from REEF3D-CFD-FWH-Box.dat,  p_dip = -1/(4 pi) rhat.F_fwh / r^2
# isolates end-cap / non-compactness errors from the REEF3D body-force bias
import math, sys, bisect
T0 = float(sys.argv[1]) if len(sys.argv)>1 else 25.0
T1 = float(sys.argv[2]) if len(sys.argv)>2 else 1e9
obs = [(0,50,0),(0,0,50),(50,0,0),(0,-50,0)]
d={}
for l in open("REEF3D_CFD_Acoustics/REEF3D-CFD-FWH-Box.dat"):
    s=l.split()
    if len(s)<13: continue
    try: v=[float(x) for x in s[:13]]
    except: continue
    d[v[0]]=v[1:]
tb=sorted(d); B=[d[t] for t in tb]
tF=[];FF=[]
for k in range(1,len(tb)-1):
    dM=[(B[k+1][9+q]-B[k-1][9+q])/(tb[k+1]-tb[k-1]) for q in range(3)]
    tF.append(tb[k]); FF.append([-(B[k][q]+B[k][3+q]+dM[q]) for q in range(3)])
def Fat(t):
    k=bisect.bisect_left(tF,t); k=min(max(k,1),len(tF)-1); w=(t-tF[k-1])/(tF[k]-tF[k-1])
    return [FF[k-1][i]*(1-w)+FF[k][i]*w for i in range(3)]
print(f"window {T0}..{T1}: FW-H vs compact dipole of the own box force (fluctuations)")
print(f"{'observer':>14} {'ratio':>6} {'corr':>6} {'rms diff/dipole':>15}")
for n,x in enumerate(obs):
    r=math.sqrt(sum(c*c for c in x)); rh=[c/r for c in x]
    A=[];P=[]
    for l in open(f"REEF3D_CFD_Acoustics/REEF3D-CFD-FWH-Observer-{n+1}.dat"):
        s=l.split()
        if len(s)<2: continue
        try: t,p=float(s[0]),float(s[1])
        except: continue
        if t<T0 or t>T1 or t>tF[-1]: continue
        F=Fat(t-r/1500.0); A.append(p); P.append(-sum(rh[i]*F[i] for i in range(3))/(4*math.pi*r*r))
    ma=sum(A)/len(A); mp=sum(P)/len(P); a=[x-ma for x in A]; b=[x-mp for x in P]
    ra=math.sqrt(sum(x*x for x in a)/len(a)); rb=math.sqrt(sum(x*x for x in b)/len(b))
    c=sum(x*y for x,y in zip(a,b))/len(a)/(ra*rb); e=math.sqrt(sum((x-y)**2 for x,y in zip(a,b))/len(a))/rb
    print(f"{str(x):>14} {ra/rb:6.3f} {c:6.3f} {e:15.3f}")
