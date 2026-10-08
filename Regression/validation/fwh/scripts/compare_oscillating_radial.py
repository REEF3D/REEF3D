# FW-H vs CFD pressure probe vs potential-flow dipole for the oscillating sphere
# x(t) = X0 (1 - cos wt):  p = rho a^3 Udot cos(theta) / (2 r^2),  Udot = X0 w^2 cos wt
import math, sys
rho, a, X0, w = 1000.0, 0.1, 0.01, 2*math.pi
obs = [(r,0.004,0.006) for r in (0.3,0.4,0.5,0.7,1.0,1.5)]+[(0.004,0.5,0.006)]
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

ref=fit(read("REEF3D_CFD_PressureProbe/REEF3D-CFD-Probe-Pressure-7.dat"))
print(f"gauge probe 7 (cos theta ~ 0): cos {ref[0]:.4f} sin {ref[1]:.4f}")
print(f"fit window {t0}..{t1} s; cos/sin coefficients [Pa]; ratios to analytic (cos)")
print(f"{'r':>5} {'analytic':>9} {'FW-H cos':>9} {'FW-H sin':>9} {'probe cos':>9} {'probe sin':>9} {'FW-H/an':>8} {'probe/an':>8}")
for n,x in enumerate(obs[:-1]):
    r=math.sqrt(sum(c*c for c in x)); ct=x[0]/r
    an=rho*a**3*X0*w*w*ct/(2*r*r)
    fc,fs,_=fit(read(f"REEF3D_CFD_Acoustics/REEF3D-CFD-FWH-Observer-{n+1}.dat"))
    pc,ps,_=fit(read(f"REEF3D_CFD_PressureProbe/REEF3D-CFD-Probe-Pressure-{n+1}.dat"))
    pc-=ref[0]; ps-=ref[1]
    print(f"{r:5.2f} {an:9.4f} {fc:9.4f} {fs:9.4f} {pc:9.4f} {ps:9.4f} {fc/an:8.3f} {pc/an:8.3f}")
