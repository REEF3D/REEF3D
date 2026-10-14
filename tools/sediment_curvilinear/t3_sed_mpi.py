# P5a T3: several ranks. (a) the dam break of T1 (bed load, S 31 1, hydrostatic): curvilinear on the Cartesian
# nodes against the Cartesian run, per rank, on 2 and 4 ranks; (b) the river of T2 with the steep mound and the
# sand slide (S 90 1) on the grid turned by 30 deg: 2 and 4 ranks against one rank (sorted values; hydrostatic,
# fixed time step)
import sys, numpy as np
from p5lib import *
steps = int(sys.argv[1]) if len(sys.argv) > 1 else 30
sed = "S 10 1\nS 11 2\nS 20 0.0005\nS 15 3\nS 17 1.0\nS 35 20.0\nS 41 1\nS 43 1\nS 44 1\nS 31 1\n"
# (a)
nx, ny, nz, dx = 32, 24, 6, 0.125
for np_ in ([2, 4] if os.environ.get("PARTA","1")=="1" else []):
    ctl = control(nx, ny, nz, nx*dx, ny*dx, 0.5, np_=np_)
    ctrl = "A 10 5\nA 520 0\nN 47 0.3\nN 48 0\nN 49 0.01\nF 60 0.2\nF 72 1.0 1.75 1.0 2.0 0.35\nM 10 %d\nP 10 0\nW 22 -9.81\n" % np_ + sed
    res = {}
    for name, g in [("cart", None), ("curv", nodes(nx, ny, dx, "id"))]:
        d = os.path.join(S5, "runs/t3a_%s_%d" % (name, np_))
        grid = None if g is None else (g[0], g[1], np.zeros_like(g[0]))
        rc = run(d, ctl, ctrl + ("G 50 2\n" if g is not None else ""), grid=grid, np_=np_, steps=steps)
        res[name] = [final_state(d, r) for r in range(np_)]
    print("(a) np %d curv vs cart: " % np_ + "  ".join("%s %.1e" % (f, max(np.abs(res["curv"][r][f]-res["cart"][r][f]).max() for r in range(np_)))
          for f in ["bed", "eta", "U", "V"]), "| bed change max %.2e" % max(np.abs(s["bed"]).max() for s in res["cart"]))
# (b)
nx, ny, nz, dx = 40, 12, 6, 0.1
W, h = ny*dx, 0.25
bc = "C 11 1\nC 12 21\nC 13 21\nC 14 2\nC 15 21\nC 16 3\n"
I, J = np.meshgrid(np.arange(nx+1), np.arange(ny+1), indexing="ij")
Z = np.maximum(0.1 - np.hypot(I*dx-1.5, J*dx-0.5*W)*0.6, 0.0)
X, Y = nodes(nx, ny, dx, "rot", 2.7e5, 7.04e6, 30.0)
ref = None
for np_ in [1, 2, 4]:
    ctl = control(nx, ny, nz, nx*dx, W, 0.5, bc=bc, np_=np_)
    ctrl = ("A 10 5\nA 519 1\nA 520 0\nA 514 0\nA 524 1\nB 50 0.01\nB 60 1\nW 10 %g\nF 60 %g\nF 62 %g\nN 47 0.3\nN 48 0\nN 49 0.01\nM 10 %d\nP 10 0\nW 22 -9.81\n"
            % (0.4*W*h, h, h, np_)) + sed + os.environ.get("SLIDE", "S 90 1\nS 81 8\nS 82 0\nS 35 50\n").replace(";", "\n")
    d = os.path.join(S5, "runs/t3b_np%d" % np_)
    rc = run(d, ctl, ctrl + "G 50 2\n", grid=(X, Y, Z), np_=np_, steps=steps)
    v = {f: np.sort(np.concatenate([final_state(d, r)[f] for r in range(np_)])) for f in ["bed", "eta", "U", "V"]}
    if ref is None:
        ref = v
        print("(b) np 1: rc %d  bed range %.4f %.4f" % (rc, v["bed"].min(), v["bed"].max()), log_lines(d, "!!!")[:2])
        continue
    print("(b) np %d vs np 1 (turned by 30 deg, sand slide): rc %d " % (np_, rc) + "  ".join("%s %.1e" % (f, np.abs(v[f]-ref[f]).max()) for f in v))
