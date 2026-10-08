# As al_reference.py, but the force field is the smeared actuator-line distribution of
# sixdof_actuator_disk::weights (Gaussian tubes of width eps around the blade lines, Hough-Ordway
# radial shape) on a fine grid, normalised to T and Q, instead of thin lines.
import math, sys
rho, n, D, Z, sense = 1000.0, 1.0, 1.0, 4, 1
T, Q = 0.2*rho*n*n*D**4, 0.03*rho*n*n*D**5
R, Rh, eps = 0.5*D, 0.1*D, 0.1
h = float(sys.argv[1]) if len(sys.argv)>1 else 0.025
axis = (-1.0,0.0,0.0); e1 = (0.0,0.0,1.0)
def cross(a,b): return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
e2 = tuple(sense*c for c in cross(axis,e1))
def shape(r):
    rh = Rh/R; s = (r/R-rh)/(1.0-rh); return s*math.sqrt(max(1.0-s,0.0))

# grid points of the rotor region (x along the axis, y, z)
pts=[]
na=int(round(3*eps/h)); nr=int(round(R/h))
for i in range(-na,na+1):
    for j in range(-nr,nr+1):
        for k in range(-nr,nr+1):
            x,y,z=i*h,j*h,k*h
            r=math.hypot(y,z)
            if Rh<r<R: pts.append((x,y,z,r))
V=h**3

def field(phase):
    # weights wa, wt and et of every point, normalised as in the coupling (fl = 1)
    F=[]; SA=0.0; ST=0.0
    for (x,y,z,r) in pts:
        a=-x                                # coordinate along axis = (-1,0,0)
        rv=(0.0,y,z)
        th=math.atan2(sum(rv[i]*e2[i] for i in range(3)), sum(rv[i]*e1[i] for i in range(3)))
        fa=0.0
        for q in range(Z):
            dth=th-phase-2*math.pi*q/Z
            if math.cos(dth)<=0: continue
            s=r*math.sin(dth); fa+=math.exp(-(a*a+s*s)/eps**2)
        sh=shape(r); wa=sh*fa; wt=sh/(r/R)*fa
        et=tuple(sense*c for c in cross(axis,(0.0,y/r,z/r)))
        F.append((x,y,z,wa,wt,et)); SA+=wa*V; ST+=wt*r*V
    return [((x,y,z),tuple((-T*wa/SA*axis[i] + Q*wt/ST*et[i])*V for i in range(3))) for (x,y,z,wa,wt,et) in F]

obs=[(0.003,1.0,0.004),(0.003,0.004,1.5),(0.003,3.0,0.004),(-1.5,0.003,0.004),(0.5,1.2,0.3)]
w=2*math.pi*Z*n
NT=16
tt=[k*(1.0/(Z*n))/NT for k in range(NT)]
P=[[0.0]*NT for _ in obs]
for k,t in enumerate(tt):
    fl=field(2*math.pi*n*t)
    for o,xo in enumerate(obs):
        p=0.0
        for (y,f) in fl:
            rv=(xo[0]-y[0],xo[1]-y[1],xo[2]-y[2]); rr=math.sqrt(rv[0]**2+rv[1]**2+rv[2]**2)
            p+=(f[0]*rv[0]+f[1]*rv[1]+f[2]*rv[2])/(4*math.pi*rr**3)
        P[o][k]=p
print(f"smeared actuator line (eps {eps}, grid {h}, {len(pts)} points): mean, BPF amplitude, phase [deg]")
for o,xo in enumerate(obs):
    m=sum(P[o])/NT
    c=2.0/NT*sum(P[o][k]*math.cos(w*tt[k]) for k in range(NT)); s=2.0/NT*sum(P[o][k]*math.sin(w*tt[k]) for k in range(NT))
    print(f"{str(xo):>22}  {m:8.3f} {math.hypot(c,s):8.3f} {math.degrees(math.atan2(-s,c)):7.1f}")
