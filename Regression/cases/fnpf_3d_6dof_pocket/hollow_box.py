# Architect: Hans Bihs
# Box with a closed inner cavity for the FNPF pocket regression case (fnpf_3d_6dof_pocket): the box
# 0.6 x 0.6 m, z 0.8 .. 1.15 (still water at 1.0), centre (5, 1); the cavity 0.11 x 0.1 x 0.08 m at
# z 0.86 .. 0.94 around the grid nodes at (5.0, 1.0) (off-centre in x: the ray there misses the diagonals of the
# cavity faces). A vertical ray through the cavity crosses the hull four times:
# its nodes are fluid that reaches neither the free surface nor the rest of the fluid (close_pockets).
import numpy as np, sys
def box(lo, hi, inward=False):
    x0, y0, z0 = lo; x1, y1, z1 = hi
    v = [np.array(p, float) for p in [(x0,y0,z0),(x1,y0,z0),(x1,y1,z0),(x0,y1,z0),(x0,y0,z1),(x1,y0,z1),(x1,y1,z1),(x0,y1,z1)]]
    q = [(0,3,2,1),(4,5,6,7),(0,1,5,4),(1,2,6,5),(2,3,7,6),(3,0,4,7)]   # outward
    t = []
    for a, b, c, d in q:
        f = [(v[a], v[b], v[c]), (v[a], v[c], v[d])]
        t += [(p[0], p[2], p[1]) for p in f] if inward else f
    return t
tri = box((4.7, 0.7, 0.8), (5.3, 1.3, 1.15)) + box((4.95, 0.95, 0.86), (5.06, 1.05, 0.94), inward=True)
with open(sys.argv[1], 'w') as f:
    f.write('solid hollow\n')
    for t in tri:
        n = np.cross(t[1] - t[0], t[2] - t[0]); n /= np.linalg.norm(n)
        f.write('  facet normal %.9e %.9e %.9e\n    outer loop\n' % tuple(n))
        for v in t: f.write('      vertex %.9e %.9e %.9e\n' % tuple(v))
        f.write('    endloop\n  endfacet\n')
    f.write('endsolid hollow\n')
