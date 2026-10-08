import math, glob, sys
rho,a,X0,w=1000.0,0.1,0.01,2*math.pi
t0,t1=0.25,1.0
def read(fn,col):
    d=[]
    for l in open(fn):
        s=l.split()
        try: d.append((float(s[0]),float(s[col])))
        except: pass
    return d
def fit(d):
    import itertools
    M=[[0.0]*3 for _ in range(3)]; b=[0.0]*3
    for t,p in d:
        if t<t0 or t>t1: continue
        f=(math.cos(w*t),math.sin(w*t),1.0)
        for i in range(3):
            b[i]+=f[i]*p
            for j in range(3): M[i][j]+=f[i]*f[j]
    def det(m): return (m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])-m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])+m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]))
    D=det(M)
    return [det([[b[r] if cc==c else M[r][cc] for cc in range(3)] for r in range(3)])/D for c in range(3)]
files=sorted(glob.glob("REEF3D_CFD_ProbePoint/*.dat"), key=lambda f:int(''.join(ch for ch in f.split('-')[-1] if ch.isdigit())))
rs=(0.15,0.2,0.3,0.5)
g=fit(read(files[4],4))
print("u_x on the x axis: sin coefficient vs potential flow X0 w a^3/r^3; p: cos coefficient (minus probe at theta=90) vs rho a^3 X0 w^2/(2 r^2)")
print(f"{'r':>5} {'u an':>9} {'u CFD':>9} {'ratio':>6} {'p an':>8} {'p CFD':>8} {'ratio':>6}")
for n,r in enumerate(rs):
    uc=fit(read(files[n],1)); pc=fit(read(files[n],4))
    ua=X0*w*a**3/r**3; pa=rho*a**3*X0*w*w/(2*r*r)
    print(f"{r:5.2f} {ua:9.5f} {uc[1]:9.5f} {uc[1]/ua:6.3f} {pa:8.4f} {pc[0]-g[0]:8.4f} {(pc[0]-g[0])/pa:6.3f}")
print("FW-H (cos coefficient) vs analytic and pressure probes P 64 (minus probe 7)")
def read2(fn):
    return read(fn,1)
g7=fit(read2("REEF3D_CFD_PressureProbe/REEF3D-CFD-Probe-Pressure-7.dat"))
for n,r in enumerate((0.3,0.4,0.5,0.7,1.0,1.5)):
    pa=rho*a**3*X0*w*w/(2*r*r)
    f=fit(read2(f"REEF3D_CFD_Acoustics/REEF3D-CFD-FWH-Observer-{n+1}.dat"))
    q=fit(read2(f"REEF3D_CFD_PressureProbe/REEF3D-CFD-Probe-Pressure-{n+1}.dat"))
    print(f"{r:5.2f}  FW-H/an {f[0]/pa:6.3f}  probe/an {(q[0]-g7[0])/pa:6.3f}  FW-H sin {f[1]:8.4f} probe sin {q[1]-g7[1]:8.4f}")
