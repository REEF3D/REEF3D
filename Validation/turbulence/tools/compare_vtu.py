# Architect: Hans Bihs
"""compare_vtu.py <serial run> <mpi run> [module]: compare the last VTU/VTP output of a serial and a
multi-rank run point by point (duplicate points on rank interfaces are averaged)."""
import sys, glob, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import r3read
def last_pieces(run, sub):
    fs = sorted(glob.glob(f'{run}/{sub}/*-0000*.vtu'))
    last = fs[-1].rsplit('-', 1)[0]
    return sorted(glob.glob(last + '-*.vtu'))
def merged(run, sub):
    acc = {}
    for f in last_pieces(run, sub):
        d = r3read.read_vtk(f); P = d['points']
        # key: (x, y, vertical index in the column) - z of sigma-grid points depends on the solution
        xy = np.round(P[:, :2], 5); key = {}
        cols = {}
        for idx, (x, y) in enumerate(map(tuple, xy)): cols.setdefault((x, y), []).append(idx)
        for (x, y), ids in cols.items():
            ids = sorted(ids, key=lambda i: P[i, 2])
            for kk, i in enumerate(ids): key[i] = (x, y, kk)
        for name, v in d.items():
            if name in ('points', 'connectivity', 'offsets', 'types') or len(v) != len(P): continue
            for i, val in enumerate(v):
                acc.setdefault(name, {}).setdefault(key[i], []).append(val)
    return {n: {p: np.mean(vs, axis=0) for p, vs in m.items()} for n, m in acc.items()}
def compare(a, b, sub='REEF3D_NHFLOW_VTU'):
    A, B = merged(a, sub), merged(b, sub); out = {}
    for n in A:
        if n not in B: continue
        keys = sorted(set(A[n]) & set(B[n]))
        x = np.array([A[n][k] for k in keys]); y = np.array([B[n][k] for k in keys])
        sc = np.max(np.abs(x)) + 1e-30
        out[n] = (np.max(np.abs(x - y)) / sc, np.mean(np.abs(x - y)) / sc, sc, len(keys))
    return out
if __name__ == '__main__':
    for n, (mx, mean, sc, npnt) in compare(sys.argv[1], sys.argv[2]).items():
        print(f"{n:12s} max|diff|/max {mx:9.2e}   mean|diff|/max {mean:9.2e}   (max {sc:.3g}, {npnt} points)")
