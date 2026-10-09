# KP505 resolved (direct forcing) open-water run: KT, KQ from the 6DOF forces (surface integral) and
# from the FW-H box momentum balance against the open-water curves, and the BPF content of FW-H and
# the pressure probes. n (X 211), U (U 50), the observers (U 30) and probes (P 64) are read from
# ctrl.txt of the run directory.
#   python3 kp505_analysis.py <t0> [t1]
import math, sys
T0 = float(sys.argv[1]) if len(sys.argv)>1 else 0.4
T1 = float(sys.argv[2]) if len(sys.argv)>2 else 1e9
rho, D, Z = 1000.0, 0.25, 5

# KP505 open-water curves, (KT, 10 KQ), read from the plot of the user (+-0.005)
EXP = {0.1:(0.440,0.641),0.2:(0.400,0.589),0.3:(0.358,0.534),0.4:(0.312,0.477),0.5:(0.266,0.417),
       0.6:(0.219,0.357),0.7:(0.169,0.293),0.8:(0.117,0.226),0.9:(0.058,0.149)}
def exp_at(J):
    js=sorted(EXP)
    if J<=js[0] or J>=js[-1]: return None
    for a,b in zip(js,js[1:]):
        if a<=J<=b:
            w=(J-a)/(b-a); return tuple(EXP[a][i]*(1-w)+EXP[b][i]*w for i in range(2))

n=U=None; obs=[]; nprobe=0
for l in open("ctrl.txt"):
    s=l.split()
    if len(s)<2: continue
    if s[0]=="X" and s[1]=="211": n=abs(float(s[2]))/(2*math.pi)
    if s[0]=="U" and s[1]=="50": U=float(s[2])
    if s[0]=="U" and s[1]=="30": obs.append(tuple(float(x) for x in s[2:5]))
    if s[0]=="P" and s[1]=="64": nprobe+=1
J = U/(n*D)

def table(fn, ncol):
    d={}
    for l in open(fn):
        s=l.split()
        if len(s)<ncol: continue
        try: v=[float(x) for x in s[:ncol]]
        except: continue
        d[v[0]]=v[1:]
    t=sorted(d); return t,[d[x] for x in t]

tf,F = table("REEF3D_CFD_6DOF/REEF3D_6DOF_forces_0.dat",13)
sel=[k for k,t in enumerate(tf) if T0<=t<=T1]
mean=lambda i: sum(F[k][i] for k in sel)/len(sel)
Fx,Fy,Fz,Mx = mean(0),mean(1),mean(2),mean(3)
KT=-Fx/(rho*n*n*D**4); KQ=Mx/(rho*n*n*D**5)
print(f"window {T0}..{min(T1,tf[-1]):.4f} s ({len(sel)} steps), n = {n:.3f} 1/s, U = {U} m/s, J = {J:.3f}")
e=exp_at(J)
ex = f"   experiment KT = {e[0]:.3f}, 10 KQ = {e[1]:.3f}" if e else ""
print(f"6DOF force (surface integral): Fx = {Fx:.3f} N, side force Fy = {Fy:.3f} Fz = {Fz:.3f} N, Mx = {Mx:.4f} Nm")
print(f"  KT = {KT:.4f}, 10 KQ = {10*KQ:.4f}, eta0 = {J*KT/(2*math.pi*KQ) if KQ else float('nan'):.3f}{ex}")

tb,B = table("REEF3D_CFD_Acoustics/REEF3D-CFD-FWH-Box.dat",13)
Fb=[[],[],[]]
for k in range(1,len(tb)-1):
    if not (T0<=tb[k]<=T1): continue
    for q in range(3):
        dM=(B[k+1][9+q]-B[k-1][9+q])/(tb[k+1]-tb[k-1])
        Fb[q].append(-(B[k][q]+B[k][3+q]-B[k][6+q]+dM))
Fbm=[sum(v)/len(v) for v in Fb]
print(f"box momentum balance: Fx = {Fbm[0]:.3f} N (KT_box = {-Fbm[0]/(rho*n*n*D**4):.4f}), side force Fy = {Fbm[1]:.3f} Fz = {Fbm[2]:.3f} N")

w=2*math.pi*Z*n
def fit(d):
    rows=[(t,p) for t,p in d if T0<=t<=T1]
    m=5; A=[[0.0]*m for _ in range(m)]; b=[0.0]*m
    for t,p in rows:
        f=[1.0,math.cos(w*t),math.sin(w*t),math.cos(2*w*t),math.sin(2*w*t)]
        for i in range(m):
            b[i]+=f[i]*p
            for j in range(m): A[i][j]+=f[i]*f[j]
    for i in range(m):
        piv=max(range(i,m),key=lambda r:abs(A[r][i])); A[i],A[piv]=A[piv],A[i]; b[i],b[piv]=b[piv],b[i]
        for r in range(m):
            if r!=i:
                fct=A[r][i]/A[i][i]
                for j in range(m): A[r][j]-=fct*A[i][j]
                b[r]-=fct*b[i]
    c=[b[i]/A[i][i] for i in range(m)]
    return c[0], math.hypot(c[1],c[2]), math.degrees(math.atan2(-c[2],c[1])), math.hypot(c[3],c[4])
def rows(fn):
    d=[]
    for l in open(fn):
        s=l.split()
        if len(s)<2: continue
        try: d.append((float(s[0]),float(s[1])))
        except: pass
    return d
print(f"BPF = {Z*n:.1f} Hz; amplitude [Pa] and phase [deg] of BPF, amplitude of 2 BPF (probe k = observer k for k <= {nprobe})")
print(f"{'observer':>22} | {'FW-H BPF':>9} {'ph':>6} {'2BPF':>7} | {'probe BPF':>9} {'ph':>6} {'2BPF':>7}")
for k,x in enumerate(obs):
    fw=fit(rows(f"REEF3D_CFD_Acoustics/REEF3D-CFD-FWH-Observer-{k+1}.dat"))
    pr=""
    if k<nprobe:
        q=fit(rows(f"REEF3D_CFD_PressureProbe/REEF3D-CFD-Probe-Pressure-{k+1}.dat"))
        pr=f"{q[1]:9.4f} {q[2]:6.1f} {q[3]:7.4f}"
    print(f"{str(x):>22} | {fw[1]:9.4f} {fw[2]:6.1f} {fw[3]:7.4f} | {pr}")
