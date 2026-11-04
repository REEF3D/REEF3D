# T9b: the curvilinear diffusion (sflow_ediff::curv_stress) against the Cartesian operator on the same nodes: the dam
# break of T1 (closed basin, box dam, radial spreading) without and with A 212 1, W 2 0.05; the effect of the
# diffusion on eta (with - without) on the curvilinear and on the Cartesian grid, and their difference (the two
# operators agree to second order, not to round-off: the transposed part and the walls differ)
# usage: python3 t9b_diffusion_cartesian.py [steps] [nu]
import sys, numpy as np
from p1lib import *

n, dx = 120, 0.25
hdam = float(os.environ.get("HDAM", "2.0"))
steps = int(sys.argv[1]) if len(sys.argv) > 1 else 160
nu = float(sys.argv[2]) if len(sys.argv) > 2 else 0.05
control = "B 1 %g\nB 10 0.0 %g 0.0 %g 0.0 0.25\n\nM 10 1\nM 20 2\n" % (dx, n*dx, n*dx)
ctrl = ("A 10 2\nA 210 3\nA 211 0\nA 220 0\nA 240 1\nF 60 1.0\nF 72 12.0 18.0 12.0 18.0 HDAM\nN 41 30.0\nN 47 0.3\n"
        "N 48 0\nN 49 0.01\nM 10 1\nW 22 -9.81\n").replace("HDAM", str(hdam))
I, J = np.meshgrid(np.arange(n+1), np.arange(n+1), indexing="ij")
grid = (I*dx, J*dx, np.zeros(I.shape))
R = {}
for g in ("cart", "curv"):
    for m, extra in (("none", ""), ("visc", "A 212 1\nW 2 %g\n" % nu)):
        c = ctrl + extra + ("G 50 2\n" if g == "curv" else "")
        d = os.path.join(S, "runs/t9b_%s_%s" % (g, m))
        rc = run(d, control, c, grid=grid if g == "curv" else None, steps=steps)
        R[(g, m)] = final_state(d)["eta"]
        print("%-4s %-4s rc=%d eta %.5f .. %.5f" % (g, m, rc, R[(g, m)].min(), R[(g, m)].max()))
dc = R[("cart", "visc")] - R[("cart", "none")]
dv = R[("curv", "visc")] - R[("curv", "none")]
print("without diffusion: curv vs cart %.2e" % np.abs(R[("curv", "none")] - R[("cart", "none")]).max())
print("effect of the diffusion: cartesian %.3e, curvilinear %.3e, difference %.3e (%.1f %%)"
      % (np.abs(dc).max(), np.abs(dv).max(), np.abs(dv-dc).max(), 100*np.abs(dv-dc).max()/np.abs(dc).max()))
