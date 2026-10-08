# Open-water actuator line: FW-H and pressure probes against the incompressible dipole field of the
# rotating line forces, p = 1/(4 pi) sum f.(x - y)/|x - y|^3 (f: force on the fluid; the linear
# pressure of a momentum source, del^2 p = div f). Harmonic fit: mean, BPF and 2 BPF.
import math, sys
T0 = float(sys.argv[1]) if len(sys.argv)>1 else 2.0
T1 = float(sys.argv[2]) if len(sys.argv)>2 else 1e9

# ship.dat: D = 1, hub 0.2, KT 0.2, KQ 0.03, n = 1, Z = 4, sense 1; body yawed by 180 deg
rho, n, D, Z, sense = 1000.0, 1.0, 1.0, 4, 1
T = 0.2*rho*n*n*D**4
Q = 0.03*rho*n*n*D**5
R, Rh = 0.5*D, 0.1*D
axis = (-1.0,0.0,0.0)                 # thrust direction on the body
e1 = (0.0,0.0,1.0)                    # blade 0 at t = 0
def cross(a,b): return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
e2 = tuple(sense*c for c in cross(axis,e1))
centre = (0.0,0.0,0.0)

# radial stations, Hough-Ordway shape as in sixdof_actuator_disk::weights
NS = 60
rs_ = [Rh + (R-Rh)*(j+0.5)/NS for j in range(NS)]
dr = (R-Rh)/NS
def shape(r):
    rh = Rh/R; s = (r/R-rh)/(1.0-rh); return s*math.sqrt(max(1.0-s,0.0))
wa = [shape(r) for r in rs_]
wt = [shape(r)/(r/R) for r in rs_]
Ca = T/(Z*sum(w*dr for w in wa))                 # axial line density: -T axis in total
Ct = Q/(Z*sum(w*r*dr for w,r in zip(wt,rs_)))    # tangential: torque sense Q about the axis

def p_ref(x,t):
    p = 0.0
    for q in range(Z):
        ph = 2.0*math.pi*n*t + 2.0*math.pi*q/Z
        d = tuple(math.cos(ph)*e1[i] + math.sin(ph)*e2[i] for i in range(3))
        et = tuple(sense*c for c in cross(axis,d))
        for r,a_,t_ in zip(rs_,wa,wt):
            y = tuple(centre[i] + r*d[i] for i in range(3))
            f = tuple((-Ca*a_*axis[i] + Ct*t_*et[i])*dr for i in range(3))
            rv = tuple(x[i]-y[i] for i in range(3)); rr = math.sqrt(sum(c*c for c in rv))
            p += sum(f[i]*rv[i] for i in range(3))/(4.0*math.pi*rr**3)
    return p

def rows(fn):
    d=[]
    for l in open(fn):
        s=l.split()
        if len(s)<2: continue
        try: d.append((float(s[0]),float(s[1])))
        except: pass
    return d

w = 2.0*math.pi*Z*n
def fit(d):
    # p = c0 + sum_k a_k cos(k w t) + b_k sin(k w t), k = 1, 2 (least squares)
    rows_=[(t,p) for t,p in d if T0<=t<=T1]
    m=5; A=[[0.0]*m for _ in range(m)]; b=[0.0]*m
    for t,p in rows_:
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
    return c[0], math.hypot(c[1],c[2]), math.atan2(-c[2],c[1]), math.hypot(c[3],c[4]), len(rows_)

obs = [(0.003,1.0,0.004),(0.003,0.004,1.5),(0.003,3.0,0.004),(-1.5,0.003,0.004),(0.5,1.2,0.3),(0.0,20.0,0.0)]
print(f"window {T0}..{T1} s, BPF = {Z*n} Hz; mean, BPF amplitude and phase [deg], 2BPF amplitude [Pa]")
print(f"{'observer':>22} | {'ref mean':>9} {'ref BPF':>8} {'ph':>6} | {'FW-H mean':>9} {'BPF':>8} {'ph':>6} {'2BPF':>7} | {'probe BPF':>9} {'ph':>6}")
for k,x in enumerate(obs):
    tt=[T0+i*(1.0/(Z*n))/40.0 for i in range(80)]
    ref=fit([(t,p_ref(x,t)) for t in tt])
    fw=fit(rows(f"REEF3D_CFD_Acoustics/REEF3D-CFD-FWH-Observer-{k+1}.dat"))
    pr=""
    if k<5:
        q=fit(rows(f"REEF3D_CFD_PressureProbe/REEF3D-CFD-Probe-Pressure-{k+1}.dat"))
        pr=f"{q[1]:9.3f} {math.degrees(q[2]):6.1f}"
    print(f"{str(x):>22} | {ref[0]:9.3f} {ref[1]:8.3f} {math.degrees(ref[2]):6.1f} | {fw[0]:9.3f} {fw[1]:8.3f} {math.degrees(fw[2]):6.1f} {fw[3]:7.3f} | {pr}")
