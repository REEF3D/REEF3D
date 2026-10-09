import sys, collections
tri=[]; cur=[]
for l in open(sys.argv[1]):
    s=l.split()
    if s and s[0]=='vertex':
        cur.append(tuple(round(float(x),9) for x in s[1:4]))
        if len(cur)==3: tri.append(tuple(cur)); cur=[]
idx={}; T=[]
for t in tri:
    T.append(tuple(idx.setdefault(v,len(idx)) for v in t))
edges=collections.Counter()
for a,b,c in T:
    for e in ((a,b),(b,c),(c,a)): edges[tuple(sorted(e))]+=1
cnt=collections.Counter(edges.values())
V=list(idx.keys())
vol=0.0
for t in tri:
    (x1,y1,z1),(x2,y2,z2),(x3,y3,z3)=t
    vol+=(x1*(y2*z3-y3*z2)-x2*(y1*z3-y3*z1)+x3*(y1*z2-y2*z1))/6.0
# directed edges: each must appear once in each direction for a consistently oriented closed surface
de=collections.Counter()
for a,b,c in T:
    for e in ((a,b),(b,c),(c,a)): de[e]+=1
bad=sum(1 for e,n in de.items() if n!=1 or de.get((e[1],e[0]),0)!=1)
xs=[v[0] for v in V]; ys=[v[1] for v in V]; zs=[v[2] for v in V]
print(f"triangles {len(T)}, vertices {len(V)}, edges by number of triangles {dict(cnt)}, inconsistently oriented edges {bad}")
print(f"signed volume {vol:.6e} m^3 (positive: normals outward), bbox x {min(xs):.4f}..{max(xs):.4f} y {min(ys):.4f}..{max(ys):.4f} z {min(zs):.4f}..{max(zs):.4f}")
