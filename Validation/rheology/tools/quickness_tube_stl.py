# closed quarter-annulus prism, the lifted cylinder of the quickness test (case 20), ASCII STL: python3 quickness_tube_stl.py floating.stl
# inner radius 54 mm: one cell (4 mm) outside the clay column (r 50 mm); at r = 50 mm the cells next to the body lost 16 % of the clay at the start
import math, sys
ri, ro, z0, z1 = 0.054, 0.066, -0.01, 0.13
a0, a1, N = math.radians(-15.0), math.radians(105.0), 48
out = sys.argv[1]
def P(r, a, z): return (r*math.cos(a), r*math.sin(a), z)
tris = []
def quad(a, b, c, d): tris.append((a, b, c)); tris.append((a, c, d))
for n in range(N):
    t0 = a0 + (a1 - a0)*n/N; t1 = a0 + (a1 - a0)*(n + 1)/N
    quad(P(ro, t0, z0), P(ro, t1, z0), P(ro, t1, z1), P(ro, t0, z1))     # outer, normal outward
    quad(P(ri, t1, z0), P(ri, t0, z0), P(ri, t0, z1), P(ri, t1, z1))     # inner, normal inward (to the axis)
    quad(P(ri, t0, z1), P(ro, t0, z1), P(ro, t1, z1), P(ri, t1, z1))     # top
    quad(P(ri, t1, z0), P(ro, t1, z0), P(ro, t0, z0), P(ri, t0, z0))     # bottom
for a, s in ((a0, -1), (a1, 1)):                                           # end caps
    A, B, C, D = P(ri, a, z0), P(ro, a, z0), P(ro, a, z1), P(ri, a, z1)
    if s > 0: quad(A, D, C, B)
    else: quad(A, B, C, D)
def nrm(a, b, c):
    u = [b[i] - a[i] for i in range(3)]; v = [c[i] - a[i] for i in range(3)]
    n = (u[1]*v[2] - u[2]*v[1], u[2]*v[0] - u[0]*v[2], u[0]*v[1] - u[1]*v[0]); l = math.sqrt(sum(x*x for x in n)) or 1.0
    return [x/l for x in n]
with open(out, "w") as f:
    f.write("solid tube\n")
    for t in tris:
        n = nrm(*t)
        f.write("  facet normal %.8f %.8f %.8f\n    outer loop\n" % tuple(n))
        for v in t: f.write("      vertex %.8f %.8f %.8f\n" % v)
        f.write("    endloop\n  endfacet\n")
    f.write("endsolid tube\n")
# centroid (z only matters for the motion file)
print("triangles", len(tris), "z centre", 0.5*(z0 + z1))
