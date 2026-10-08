# FW-H vs CFD pressure probe vs potential-flow dipole for the oscillating sphere
# x(t) = X0 (1 - cos wt):  p = rho a^3 Udot cos(theta) / (2 r^2),  Udot = X0 w^2 cos wt
import math, sys
rho, a, X0, w = 1000.0, 0.1, 0.01, 2*math.pi
obs = [(0.5,0,0),(0.353553,0.353553,0),(0,0.5,0),(-0.4,0,0.3),(5.0,0,0)]
t0, t1 = float(sys.argv[1]) if len(sys.argv)>1 else 0.5, float(sys.argv[2]) if len(sys.argv)>2 else 2.0

def read(fn):
    d=[]
    for l in open(fn):
        s=l.split()
        try: d.append((float(s[0]),float(s[1])))
        except: pass
    return d

def solve3(M,b):
    # Cramer
    def det(m): return (m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])-m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])+m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]))
    D=det(M); x=[]
    for c in range(3):
        Mc=[[b[r] if cc==c else M[r][cc] for cc in range(3)] for r in range(3)]
        x.append(det(Mc)/D)
    return x

def fit(d):
    M=[[0.0]*3 for _ in range(3)]; b=[0.0]*3; n=0
    for t,p in d:
        if t<t0 or t>t1: continue
        f=(math.cos(w*t),math.sin(w*t),1.0); n+=1
        for i in range(3):
            b[i]+=f[i]*p
            for j in range(3): M[i][j]+=f[i]*f[j]
    c=solve3(M,b)
    return c[0],c[1],n

# probes relative to probe 3 (cos theta = 0, p = 0): removes the floating pressure level of the closed domain
ref=fit(read("REEF3D_CFD_PressureProbe/REEF3D-CFD-Probe-Pressure-3.dat"))
print(f"probe 3 (gauge) fit: cos {ref[0]:.4f} sin {ref[1]:.4f}")
print(f"fit window {t0}..{t1} s, coefficients of cos(wt), sin(wt) [Pa]")
print(f"{'observer':>24} {'analytic':>9} {'FW-H cos':>9} {'FW-H sin':>9} {'probe cos':>9} {'probe sin':>9}")
for n,x in enumerate(obs):
    r=math.sqrt(sum(c*c for c in x)); ct=x[0]/r
    an=rho*a**3*X0*w*w*ct/(2*r*r)
    f=read(f"REEF3D_CFD_Acoustics/REEF3D-CFD-FWH-Observer-{n+1}.dat")
    fc,fs,nf=fit(f)
    pr=""
    if n<4:
        pd=read(f"REEF3D_CFD_PressureProbe/REEF3D-CFD-Probe-Pressure-{n+1}.dat")
        pc,ps,_=fit(pd); pr=f"{pc-ref[0]:9.4f} {ps-ref[1]:9.4f}"
    print(f"{str(x):>24} {an:9.4f} {fc:9.4f} {fs:9.4f} {pr}")
