# compare two AMR runs cell by cell (level, I, J): eta, U, V
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../sediment_curvilinear"))
from p5lib import amr_grids2
def tab(d, np_):
    t = {}
    for g in amr_grids2(d, np_):
        for a in range(g["nx"]):
            for b in range(g["ny"]):
                t[(g["lev"], g["blk"], g["I0"]+a, g["J0"]+b)] = (g["eta"][a, b], g["U"][a, b], g["V"][a, b])
    return t
A, B = tab(sys.argv[1], int(sys.argv[2])), tab(sys.argv[3], int(sys.argv[4]))
for lv in sorted({k[0] for k in A}):
    com = [k for k in A if k[0] == lv and k in B]
    e = [max(abs(A[k][m]-B[k][m]) for k in com) for m in range(3)]
    print("L%d %d cells: eta %.1e U %.1e V %.1e" % (lv, len(com), *e))
