# Architect: Hans Bihs
"""nhf_compare.py <run A> <run B> [x0 x1]: NHFLOW, last VTU of two runs on the same grid. Prints the largest and
the mean difference of u, w, the vertex z (free surface) and nu_t (if written) over x0 < x < x1 (default the
whole domain), relative to the maximum of run A (z: relative to the water depth range of run A)."""
import sys, os, glob, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import r3read

def last(run):
    return r3read.read_vtk(sorted(glob.glob(run + '/REEF3D_NHFLOW_VTU/*-000001.vtu'))[-1])

if __name__ == '__main__':
    a, b = last(sys.argv[1]), last(sys.argv[2])
    x = a['points'][:, 0]
    x0, x1 = (float(sys.argv[3]), float(sys.argv[4])) if len(sys.argv) > 4 else (x.min() - 1, x.max() + 1)
    m = (x > x0) & (x < x1)
    fields = [('u', a['velocity'][:, 0], b['velocity'][:, 0]), ('w', a['velocity'][:, 2], b['velocity'][:, 2])]
    if 'eddyv' in a and 'eddyv' in b:
        fields.append(('nu_t', a['eddyv'], b['eddyv']))
    for name, fa, fb in fields:
        d = np.abs(fa[m] - fb[m]); ref = np.abs(fa[m]).max()
        print(f"   {name:5s} max |A-B|/max|A| = {d.max()/ref:.3e}   mean = {d.mean()/ref:.3e}")
    za, zb = a['points'][m, 2], b['points'][m, 2]
    d = np.abs(za - zb); ref = za.max() - za.min()
    print(f"   {'z':5s} max |A-B|/depth  = {d.max()/ref:.3e}   mean = {d.mean()/ref:.3e}")
