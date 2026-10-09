# KP505 resolved (direct forcing) open-water run: KT, KQ from the 6DOF forces and from the FW-H box
# momentum balance, and the BPF content of FW-H and the pressure probes.
import math, sys, bisect
T0 = float(sys.argv[1]) if len(sys.argv)>1 else 0.4
T1 = float(sys.argv[2]) if len(sys.argv)>2 else 1e9
rho, n, D, Z, U = 1000.0, 5.0, 0.25, 5, 0.75
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
Fx=sum(F[k][0] for k in sel)/len(sel); Mx=sum(F[k][3] for k in sel)/len(sel)
Fxp=sum(F[k][6] for k in sel)/len(sel); Fxv=sum(F[k][9] for k in sel)/len(sel)
KT=-Fx/(rho*n*n*D**4); KQ=Mx/(rho*n*n*D**5)
print(f"window {T0}..{min(T1,tf[-1]):.3f} s ({len(sel)} steps), J = {J:.3f}")
print(f"6DOF force: Fx = {Fx:.3f} N (pressure {Fxp:.3f}, viscous {Fxv:.3f}), Mx = {Mx:.4f} Nm  ->  KT = {KT:.4f}, 10 KQ = {10*KQ:.4f}, eta0 = {J*KT/(2*math.pi*KQ) if KQ!=0 else float('nan'):.3f}")

tb,B = table("REEF3D_CFD_Acoustics/REEF3D-CFD-FWH-Box.dat",13)
Fb=[]
for k in range(1,len(tb)-1):
    if not (T0<=tb[k]<=T1): continue
    dM=(B[k+1][9]-B[k-1][9])/(tb[k+1]-tb[k-1])
    Fb.append(-(B[k][0]+B[k][3]-B[k][6]+dM))
Fbx=sum(Fb)/len(Fb)
print(f"box momentum balance: Fx = {Fbx:.3f} N  ->  KT_box = {-Fbx/(rho*n*n*D**4):.4f}  (ratio to 6DOF {Fbx/Fx:.3f})")

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
obs=[(0.003,0.25,0.004),(0.003,0.004,0.375),(0.1,0.2,0.05),(-0.3,0.003,0.004),(0.0,5.0,0.0)]
print(f"BPF = {Z*n} Hz; amplitude [Pa] and phase [deg] of BPF, amplitude of 2 BPF")
print(f"{'observer':>22} | {'FW-H BPF':>9} {'ph':>6} {'2BPF':>7} | {'probe BPF':>9} {'ph':>6} {'2BPF':>7}")
for k,x in enumerate(obs):
    fw=fit(rows(f"REEF3D_CFD_Acoustics/REEF3D-CFD-FWH-Observer-{k+1}.dat"))
    pr=""
    if k<4:
        q=fit(rows(f"REEF3D_CFD_PressureProbe/REEF3D-CFD-Probe-Pressure-{k+1}.dat"))
        pr=f"{q[1]:9.4f} {q[2]:6.1f} {q[3]:7.4f}"
    print(f"{str(x):>22} | {fw[1]:9.4f} {fw[2]:6.1f} {fw[3]:7.4f} | {pr}")
