# FW-H observer planes (z = 0 grids of observers in ctrl.txt, e.g. made by
# cases/15_sphere_re300/make_case.py) to a time series of legacy VTK structured-points files for
# ParaView (open fwh_<spacing>_..vtk as a group and play):
#   python3 fwh_planes_vtk.py [t0 t1 [dt]]      in the run directory
# Without arguments (e.g. Run in Spyder): the last 20 s of the output, dt = 0.33 s.
# One file per plane spacing L and frame, REEF3D_FWH_Planes/fwh_L<L>_<frame>.vtk, with the point data p' [Pa]
# (and p' r^2 or p' r, normalised by the decay of the near or far field, to show the pattern).
# The frames of an earlier call in that folder (fwh_L*.vtk) are removed first.
import math, os, sys, bisect, collections, glob
OUT = "REEF3D_FWH_Planes"
args = [float(x) for x in sys.argv[1:4]]

obs=[]
for l in open("ctrl.txt"):
    s=l.split()
    if len(s)>=5 and s[0]=="U" and s[1]=="30":
        obs.append(tuple(float(x) for x in s[2:5]))

# planes: observers on a regular x-y grid with spacing L (centre left out)
planes=collections.defaultdict(dict)
for n,(x,y,z) in enumerate(obs):
    for L in (5.0, 5000.0):
        i,j=round(x/L),round(y/L)
        if abs(x-i*L)<1e-6*L and abs(y-j*L)<1e-6*L and abs(i)<=8 and abs(j)<=8 and abs(z)<0.1 and max(abs(i),abs(j))>0:
            planes[L][(i,j)]=n
planes={L:v for L,v in planes.items() if len(v)>=200}

def signal(n):
    d={}
    for l in open(f"REEF3D_CFD_Acoustics/REEF3D-CFD-FWH-Observer-{n+1}.dat"):
        s=l.split()
        if len(s)!=2: continue
        try: d[float(s[0])]=float(s[1])
        except: pass
    t=sorted(d); return t,[d[x] for x in t]
def at(ts,ps,t):
    k=bisect.bisect_left(ts,t)
    if k<=0 or k>=len(ts): return float('nan')
    w=(t-ts[k-1])/(ts[k]-ts[k-1]); return ps[k-1]*(1-w)+ps[k]*w

if not planes:
    sys.exit("no observer planes (regular z = 0 grids of U 30) in ctrl.txt")

# frame times: from the arguments, or the last 20 s that all plane observers have
sigs={n:signal(n) for L in planes for n in planes[L].values()}
tlast=min(sigs[n][0][-1] for n in sigs)
t0 = args[0] if len(args)>=2 else tlast-20.0
t1 = min(args[1], tlast) if len(args)>=2 else tlast
dt = args[2] if len(args)>=3 else 0.33
print(f"frames {t0:.2f}..{t1:.2f} s every {dt} s (output up to t = {tlast:.2f} s)")
nframes=int(math.floor((t1-t0)/dt+1e-9))+1
os.makedirs(OUT, exist_ok=True)
for old in glob.glob(os.path.join(OUT,"fwh_L*.vtk")):
    os.remove(old)
for L,pts in planes.items():
    sig={key:sigs[n] for key,n in pts.items()}
    # remove the mean of each observer (steady part) for the animation
    mean={key:sum(p for t,p in zip(*sig[key]) if t0<=t<=t1)/max(1,sum(1 for t in sig[key][0] if t0<=t<=t1)) for key in sig}
    far = L>100.0
    for f in range(nframes):
        t=t0+f*dt
        vals=[];norm=[]
        for j in range(-8,9):
            for i in range(-8,9):
                if (i,j) in sig:
                    p=at(*sig[(i,j)],t)-mean[(i,j)]; r=math.hypot(i*L,j*L)
                    vals.append(p); norm.append(p*(r if far else r*r))
                else:
                    vals.append(0.0); norm.append(0.0)
        name=os.path.join(OUT,f"fwh_L{int(L)}_{f:04d}.vtk")
        with open(name,"w") as o:
            o.write(f"# vtk DataFile Version 3.0\nFW-H p' t = {t:.4f} s\nASCII\nDATASET STRUCTURED_POINTS\n")
            o.write(f"DIMENSIONS 17 17 1\nORIGIN {-8*L} {-8*L} 0\nSPACING {L} {L} 1\n")
            o.write(f"POINT_DATA {len(vals)}\nSCALARS p float 1\nLOOKUP_TABLE default\n")
            o.write("\n".join(f"{v:.6e}" for v in vals)+"\n")
            o.write(f"SCALARS p_{'r' if far else 'r2'} float 1\nLOOKUP_TABLE default\n")
            o.write("\n".join(f"{v:.6e}" for v in norm)+"\n")
    print(f"plane L = {L}: {len(pts)} observers, {nframes} frames {os.path.abspath(OUT)}/fwh_L{int(L)}_*.vtk")
