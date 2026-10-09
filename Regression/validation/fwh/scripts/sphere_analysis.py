# Sphere in a uniform flow (vortex shedding): drag, lift and Strouhal number, and FW-H against the
# Curle compact dipole of the body force, from the near field (1/r^2) to the acoustic far field (1/r).
#   python3 sphere_analysis.py [t0] [t1]        in the run directory (default: the second half of the run)
# Observers (U 30), mean flow (U 50), viscosity (W 2), speed of sound (U 11) are read from ctrl.txt.
# Curle: p(x,t) = -1/(4 pi) [ rhat.F/r^2 + rhat.dF/dt/(c r) ]_(t - r/c), F the force of the fluid on
# the sphere, either from the 6DOF output (surface integral) or from the FW-H box momentum balance
# F = -(Lp + Lm + dM/dt) (REEF3D-CFD-FWH-Box.dat, the dipole the FW-H integral radiates).
import math, sys, bisect
T0 = float(sys.argv[1]) if len(sys.argv)>1 else None
T1 = float(sys.argv[2]) if len(sys.argv)>2 else 1e9
rho, D = 1000.0, 1.0
c0, U, nu = 1500.0, None, None
obs=[]
for l in open("ctrl.txt"):
    s=l.split()
    if len(s)<3: continue
    if s[0]=="U" and s[1]=="11": c0=float(s[2])
    if s[0]=="U" and s[1]=="50": U=float(s[2])
    if s[0]=="W" and s[1]=="2": nu=float(s[2])
    if s[0]=="U" and s[1]=="30": obs.append(tuple(float(x) for x in s[2:5]))
A = math.pi*D*D/4.0
q = 0.5*rho*U*U*A

def table(fn, ncol):
    d={}
    for l in open(fn):
        s=l.split()
        if len(s)<ncol: continue
        try: v=[float(x) for x in s[:ncol]]
        except: continue
        d[v[0]]=v[1:]
    t=sorted(d); return t,[d[x] for x in t]

def rows(fn):
    d={}
    for l in open(fn):
        s=l.split()
        if len(s)<2: continue
        try: d[float(s[0])]=float(s[1])
        except: pass
    return sorted(d.items())

# forces: 6DOF surface integral and box momentum balance
tf,F6 = table("REEF3D_CFD_6DOF/REEF3D_6DOF_forces_0.dat",13)
F6=[v[0:3] for v in F6]
tb,B = table("REEF3D_CFD_Acoustics/REEF3D-CFD-FWH-Box.dat",13)
tB=[]; FB=[]
for k in range(1,len(tb)-1):
    dM=[(B[k+1][9+i]-B[k-1][9+i])/(tb[k+1]-tb[k-1]) for i in range(3)]
    tB.append(tb[k]); FB.append([-(B[k][i]+B[k][3+i]+dM[i]) for i in range(3)])

tend = min(tf[-1], tB[-1])
if T0 is None:
    T0 = 0.5*tend
nwin = sum(1 for x in tf if T0<=x<=T1)
if nwin < 50:
    sys.exit(f"window {T0}..{min(T1,tend)} s has {nwin} samples: the output reaches t = {tend:.3f} s "
             f"(6DOF {tf[-1]:.3f} s, Box.dat {tB[-1]:.3f} s); give a window inside the run, e.g. "
             f"python3 sphere_analysis.py {0.5*tend:.1f} {tend:.1f}")

def window(t,V):
    return [(x,v) for x,v in zip(t,V) if T0<=x<=T1]

def interp(t,V,x):
    k=bisect.bisect_left(t,x); k=min(max(k,1),len(t)-1)
    w=(x-t[k-1])/(t[k]-t[k-1]); return [V[k-1][i]*(1-w)+V[k][i]*w for i in range(len(V[0]))]

def dft_peak(t,v,fmin=0.02,fmax=1.0,nf=2000):
    # peak of the amplitude spectrum of a non-uniformly sampled signal (trapezoidal DFT)
    m=sum(v)/len(v); v=[x-m for x in v]; T=t[-1]-t[0]; best=(0,0)
    for i in range(nf):
        f=fmin+(fmax-fmin)*i/(nf-1); w=2*math.pi*f; c=s=0.0
        for k in range(1,len(t)):
            dt=t[k]-t[k-1]
            c+=0.5*(v[k]*math.cos(w*t[k])+v[k-1]*math.cos(w*t[k-1]))*dt
            s+=0.5*(v[k]*math.sin(w*t[k])+v[k-1]*math.sin(w*t[k-1]))*dt
        a=2*math.hypot(c,s)/T
        if a>best[1]: best=(f,a)
    return best

def harmonic(t,v,f):
    # least-squares amplitude and phase of the component at f (with a mean)
    w=2*math.pi*f; M=[[0.0]*3 for _ in range(3)]; b=[0.0]*3
    for x,y in zip(t,v):
        g=(1.0,math.cos(w*x),math.sin(w*x))
        for i in range(3):
            b[i]+=g[i]*y
            for j in range(3): M[i][j]+=g[i]*g[j]
    for i in range(3):
        p=max(range(i,3),key=lambda r:abs(M[r][i])); M[i],M[p]=M[p],M[i]; b[i],b[p]=b[p],b[i]
        for r in range(3):
            if r!=i:
                fc=M[r][i]/M[i][i]
                for j in range(3): M[r][j]-=fc*M[i][j]
                b[r]-=fc*b[i]
    c=[b[i]/M[i][i] for i in range(3)]
    return math.hypot(c[1],c[2]), math.atan2(-c[2],c[1])

print(f"Re = {U*D/nu:.0f}, window {T0}..{min(T1,tf[-1]):.1f} s")
for name,(t,V) in (("6DOF surface integral",(tf,F6)),("box momentum balance",(tB,FB))):
    W=window(t,V); ts=[x for x,_ in W]
    m=[sum(v[i] for _,v in W)/len(W) for i in range(3)]
    ang=math.atan2(m[2],m[1])
    lat=[v[1]*math.cos(ang)+v[2]*math.sin(ang) for _,v in W]
    f,a=dft_peak(ts,lat)
    print(f"{name:>22}: Cd = {m[0]/q:.4f}, mean Cl = {math.hypot(m[1],m[2])/q:.4f} at {math.degrees(ang):.0f} deg (y -> z),"
          f" Cl' amplitude {a/q:.4f}, St = {f*D/U:.4f}")
    if name.startswith("box"): fst=f; ang_box=ang
print("  (Johnson & Patel 1999, Re = 300: Cd 0.656, Cl 0.069, St 0.137)")

# FW-H against Curle at the shedding frequency
k=2*math.pi*fst/c0
print(f"\nshedding frequency {fst:.4f} Hz, wavelength {c0/fst:.0f} m, near/far transition r = lambda/(2 pi) = {1/k:.0f} m")
print(f"{'observer':>28} {'r':>9} | {'|p| FW-H':>10} {'Curle box':>10} {'Curle 6DOF':>10} | {'FWH/box':>7} {'dphase':>7}")
def curle(t,V,x,tt):
    r=math.sqrt(sum(c*c for c in x)); rh=[c/r for c in x]; out=[]
    for ti in tt:
        tau=ti-r/c0
        Fv=interp(t,V,tau); Fp=interp(t,V,tau+1e-3); Fm=interp(t,V,tau-1e-3)
        dF=[(Fp[i]-Fm[i])/2e-3 for i in range(3)]
        out.append(-(sum(rh[i]*Fv[i] for i in range(3))/(r*r)+sum(rh[i]*dF[i] for i in range(3))/(c0*r))/(4*math.pi))
    return out
for n,x in enumerate(obs):
    d=[(t,p) for t,p in rows(f"REEF3D_CFD_Acoustics/REEF3D-CFD-FWH-Observer-{n+1}.dat") if T0<=t<=T1]
    if len(d)<20: print(f"{str(x):>28}  (not enough complete samples in the window)"); continue
    tt=[t for t,_ in d][::5]; pf=[p for _,p in d][::5]
    a_f,ph_f=harmonic(tt,pf,fst)
    a_b,ph_b=harmonic(tt,curle(tB,FB,x,tt),fst)
    a_6,_=harmonic(tt,curle(tf,F6,x,tt),fst)
    r=math.sqrt(sum(c*c for c in x)); dph=math.degrees((ph_f-ph_b+math.pi)%(2*math.pi)-math.pi)
    print(f"{str(tuple(round(c,2) for c in x)):>28} {r:9.3g} | {a_f:10.3e} {a_b:10.3e} {a_6:10.3e} | {a_f/a_b if a_b else float('nan'):7.3f} {dph:7.1f}")
