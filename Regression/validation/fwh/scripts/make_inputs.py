# Input files of the FW-H validation runs, written into the current directory.
#
#   python3 make_inputs.py oscillating   floating.stl (icosphere, radius 0.1 m, 5120 triangles)
#                                        6DOF_motion.dat (surge x = X0 (1 - cos wt), X0 = 0.01 m, f = 1 Hz)
#   python3 make_inputs.py hub           floating.stl (icosphere, radius 0.1 m, 1280 triangles; 11-12)
#   python3 make_inputs.py shedding [A]  6DOF_motion.dat (one lateral kick y = A (1 - cos(2 pi t/4))/2,
#                                        A = 0.1 m (default), 0 <= t <= 4 s, then at rest); the sphere of
#                                        diameter 1 is Regression/cases/cfd_3d_heave_sphere_6dof/floating.stl
import math, sys

def icosphere(a, levels, name):
    t = (1+5**0.5)/2
    V = [(-1,t,0),(1,t,0),(-1,-t,0),(1,-t,0),(0,-1,t),(0,1,t),(0,-1,-t),(0,1,-t),(t,0,-1),(t,0,1),(-t,0,-1),(-t,0,1)]
    V = [tuple(c/math.sqrt(sum(x*x for x in v)) for c in v) for v in V]
    F = [(0,11,5),(0,5,1),(0,1,7),(0,7,10),(0,10,11),(1,5,9),(5,11,4),(11,10,2),(10,7,6),(7,1,8),
         (3,9,4),(3,4,2),(3,2,6),(3,6,8),(3,8,9),(4,9,5),(2,4,11),(6,2,10),(8,6,7),(9,8,1)]
    for _ in range(levels):
        cache = {}
        def mid(i, j):
            k = (min(i,j), max(i,j))
            if k not in cache:
                m = [(V[i][c]+V[j][c])/2 for c in range(3)]
                n = math.sqrt(sum(x*x for x in m))
                V.append(tuple(x/n for x in m))
                cache[k] = len(V)-1
            return cache[k]
        NF = []
        for (i,j,k) in F:
            ab, bc, ca = mid(i,j), mid(j,k), mid(k,i)
            NF += [(i,ab,ca),(j,bc,ab),(k,ca,bc),(ab,bc,ca)]
        F = NF
    with open(name, 'w') as f:
        f.write('solid sphere\n')
        for (i,j,k) in F:
            p = [tuple(a*c for c in V[q]) for q in (i,j,k)]
            u = [p[1][c]-p[0][c] for c in range(3)]
            w = [p[2][c]-p[0][c] for c in range(3)]
            n = (u[1]*w[2]-u[2]*w[1], u[2]*w[0]-u[0]*w[2], u[0]*w[1]-u[1]*w[0])
            L = math.sqrt(sum(x*x for x in n))
            cen = [sum(pp[c] for pp in p)/3 for c in range(3)]
            if sum(n[c]*cen[c] for c in range(3)) < 0:
                p = [p[0],p[2],p[1]]
                n = tuple(-x for x in n)
            f.write(f'  facet normal {n[0]/L:.6f} {n[1]/L:.6f} {n[2]/L:.6f}\n    outer loop\n')
            for pp in p:
                f.write(f'      vertex {pp[0]:.8f} {pp[1]:.8f} {pp[2]:.8f}\n')
            f.write('    endloop\n  endfacet\n')
        f.write('endsolid sphere\n')

mode = sys.argv[1] if len(sys.argv) > 1 else ""

if mode == "oscillating":
    icosphere(0.1, 4, 'floating.stl')
    X0, w = 0.01, 2*math.pi*1.0
    with open('6DOF_motion.dat', 'w') as f:
        for i in range(int(2.2/0.0005)+1):
            t = i*0.0005
            f.write(f"{t:.6f} {X0*(1-math.cos(w*t)):.10e} 0.0\n")
elif mode == "hub":
    icosphere(0.1, 3, 'floating.stl')
elif mode == "shedding":
    A = float(sys.argv[2]) if len(sys.argv) > 2 else 0.1
    with open('6DOF_motion.dat', 'w') as f:
        for i in range(4001):
            t = i*0.001
            f.write(f"{t:.4f} 0.0 {A*0.5*(1-math.cos(2*math.pi*t/4.0)):.10e}\n")
else:
    print(__doc__ if __doc__ else "usage: make_inputs.py oscillating|shedding")
    sys.exit(1)
