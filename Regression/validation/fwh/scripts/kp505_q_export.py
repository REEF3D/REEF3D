# Q-criterion of a KP505 run from its LAGOON store, for a 3D animation in ParaView:
#   python kp505_q_export.py [t0 t1] [--out DIR]          in the run directory, then
#   pvbatch kp505_q_render.py DIR [--q 2000]
# For each output in [t0, t1] (default 0.4-0.8 s): DIR/q_####.vtk (legacy binary rectilinear grid
# of x in [-0.08, 0.45], |y|, |z| < 0.18 with Q and the axial velocity u) and DIR/body_####.vtk
# (the propeller surface at that time, from the body of the store).
import argparse, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lagoon_read import lagoon

ap = argparse.ArgumentParser()
ap.add_argument("t", nargs="*", type=float)
ap.add_argument("--out", default="qframes")
args = ap.parse_args()
t0, t1 = args.t[:2] if len(args.t) >= 2 else (0.4, 0.8)
os.makedirs(args.out, exist_ok=True)

S = lagoon("REEF3D_CFD.lagoon")
x, y, z = S.coords()
ix = np.where((x >= -0.08) & (x <= 0.45))[0]
iy = np.where(np.abs(y) <= 0.18)[0]
iz = np.where(np.abs(z) <= 0.18)[0]
xs, ys, zs = x[ix], y[iy], z[iz]
frames = [k for k, t in enumerate(S.times) if t0 - 1e-9 <= t <= t1 + 1e-9]


def write_grid(fn, Q, u):
    with open(fn, "wb") as f:
        f.write(b"# vtk DataFile Version 3.0\nKP505 Q\nBINARY\nDATASET RECTILINEAR_GRID\n")
        f.write(f"DIMENSIONS {len(xs)} {len(ys)} {len(zs)}\n".encode())
        for name, c in (("X", xs), ("Y", ys), ("Z", zs)):
            f.write(f"{name}_COORDINATES {len(c)} float\n".encode())
            f.write(c.astype(">f4").tobytes() + b"\n")
        f.write(f"POINT_DATA {Q.size}\n".encode())
        for name, a in (("Q", Q), ("u", u)):
            f.write(f"SCALARS {name} float 1\nLOOKUP_TABLE default\n".encode())
            f.write(np.ascontiguousarray(a).astype(">f4").tobytes() + b"\n")


def write_body(fn, V, T):
    with open(fn, "wb") as f:
        f.write(b"# vtk DataFile Version 3.0\nKP505 body\nBINARY\nDATASET POLYDATA\n")
        f.write(f"POINTS {len(V)} float\n".encode())
        f.write(V.astype(">f4").tobytes() + b"\n")
        cells = np.column_stack([np.full(len(T), 3), T]).astype(">i4")
        f.write(f"POLYGONS {len(T)} {cells.size}\n".encode())
        f.write(cells.tobytes() + b"\n")


for j, k in enumerate(frames):
    u = S.volume("velocity", k)[np.ix_(iz, iy, ix)]          # (z, y, x, 3)
    # velocity gradient du_c/dx_d, d = x, y, z
    G = [[np.gradient(u[..., c], xs, axis=2), np.gradient(u[..., c], ys, axis=1), np.gradient(u[..., c], zs, axis=0)]
         for c in range(3)]
    Q = np.zeros(u.shape[:3])
    for a in range(3):
        for b in range(3):
            Q -= 0.5*G[a][b]*G[b][a]
    write_grid(os.path.join(args.out, f"q_{j:04d}.vtk"), Q, u[..., 0])
    V, T = S.body(k)
    write_body(os.path.join(args.out, f"body_{j:04d}.vtk"), V, T)
    with open(os.path.join(args.out, "times.txt"), "a" if j else "w") as f:
        f.write(f"{j} {S.times[k]:.6f}\n")
    print(f"frame {j+1}/{len(frames)} t = {S.times[k]:.4f}", flush=True)
