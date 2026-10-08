# FW-H with closed vs open/averaged end cap on the same CFD run: p_variant vs p_closed and vs Curle
import math, sys, bisect
closed, variant = sys.argv[1], sys.argv[2]
T0 = float(sys.argv[3]) if len(sys.argv)>3 else 25.0
T1 = float(sys.argv[4]) if len(sys.argv)>4 else 1e9
obs = [(0,50,0),(0,0,50),(50,0,0),(0,-50,0),(0,3,0),(0,0,3)]
def rows(fn):
    d=[]
    for l in open(fn):
        s=l.split()
        if len(s)<2: continue
        try: d.append((float(s[0]),float(s[1])))
        except: pass
    return d
def stats(a,b):
    ma=sum(a)/len(a); mb=sum(b)/len(b); a=[x-ma for x in a]; b=[x-mb for x in b]
    ra=math.sqrt(sum(x*x for x in a)/len(a)); rb=math.sqrt(sum(x*x for x in b)/len(b))
    return ra/rb, sum(x*y for x,y in zip(a,b))/len(a)/(ra*rb), math.sqrt(sum((x-y)**2 for x,y in zip(a,b))/len(a))/rb
print(f"window {T0}..{T1}: variant vs closed box (fluctuations)")
print(f"{'observer':>14} {'rms ratio':>9} {'corr':>6} {'rms diff/closed':>15}")
for n,x in enumerate(obs):
    A=dict(rows(f"{closed}/REEF3D_CFD_Acoustics/REEF3D-CFD-FWH-Observer-{n+1}.dat"))
    B=rows(f"{variant}/REEF3D_CFD_Acoustics/REEF3D-CFD-FWH-Observer-{n+1}.dat")
    a=[];b=[]
    for t,p in B:
        if t<T0 or t>T1 or t not in A: continue
        a.append(p); b.append(A[t])
    if len(a)<10: print(x,"too few"); continue
    r,c,d=stats(a,b)
    print(f"{str(x):>14} {r:9.3f} {c:6.3f} {d:15.3f}")
