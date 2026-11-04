# Architect: Hans Bihs
# Floating frustum with tumblehome for the X 27 regression case: radius 0.25 m at the bottom
# (z = 0.85), 0.17 m at the still-water level (z = 1.0), vertical above to z = 1.15; centre (7, 1).
import numpy as np, sys
xc, yc, nt = 7.0, 1.0, 48
zs = [0.85, 1.0, 1.15]; rs = [0.25, 0.17, 0.17]
th = np.linspace(0, 2*np.pi, nt + 1)[:-1]
ring = lambda z, r: [np.array([xc + r*np.cos(t), yc + r*np.sin(t), z]) for t in th]
R = [ring(z, r) for z, r in zip(zs, rs)]
tri = []
for a in range(len(R) - 1):
    for k in range(nt):
        p0, p1, q0, q1 = R[a][k], R[a][(k+1) % nt], R[a+1][k], R[a+1][(k+1) % nt]
        tri += [(p0, p1, q1), (p0, q1, q0)]
cb, ct = np.array([xc, yc, zs[0]]), np.array([xc, yc, zs[-1]])
for k in range(nt):
    tri += [(cb, R[0][(k+1) % nt], R[0][k]), (ct, R[-1][k], R[-1][(k+1) % nt])]
with open(sys.argv[1], 'w') as f:
    f.write('solid frustum\n')
    for t in tri:
        n = np.cross(t[1] - t[0], t[2] - t[0]); n /= np.linalg.norm(n)
        f.write('  facet normal %.9e %.9e %.9e\n    outer loop\n' % tuple(n))
        for v in t: f.write('      vertex %.9e %.9e %.9e\n' % tuple(v))
        f.write('    endloop\n  endfacet\n')
    f.write('endsolid frustum\n')
