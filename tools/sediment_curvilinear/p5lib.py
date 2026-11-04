# helpers for the curvilinear sediment tests (G 50 with S 10, phase P5a): runs go to ./runs
# binaries: environment REEF3D_BIN, DIVEMESH_BIN (as ../sflow_curvilinear/p1lib.py)
import os, sys, re, shutil, subprocess, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../sflow_curvilinear"))
from p1lib import *
import p1lib
S5 = os.path.dirname(os.path.abspath(__file__))

def control(nx, ny, nz, L, W, H, bc="C 11 21\nC 12 21\nC 13 21\nC 14 21\nC 15 21\nC 16 3\n", np_=1):
    return "%s\nB 2 %d %d %d\nB 10 0.0 %g 0.0 %g 0.0 %g\n\nM 10 %d\nM 20 2\n" % (bc, nx, ny, nz, L, W, H, np_)

def nodes(nx, ny, dx, rot="id", X0=0.0, Y0=0.0, th=30.0, wavy=0.0):
    """node coordinates (nx+1, ny+1): id (on the Cartesian nodes), rot90 (turned by 90 deg about (X0,Y0)),
    rot (turned by th deg, shifted to (X0,Y0)), wavy (sine perturbation of amplitude wavy*dx)"""
    I, J = np.meshgrid(np.arange(nx+1), np.arange(ny+1), indexing="ij")
    s, n = I*dx, J*dx
    if wavy:
        s = s + wavy*dx*np.sin(2*np.pi*J/ny)*np.sin(np.pi*I/nx*2)
        n = n + wavy*dx*np.sin(2*np.pi*I/nx)*np.sin(np.pi*J/ny*2)
    if rot == "id":
        return s, n
    if rot == "rot90":
        return X0 - n, Y0 + s
    t = np.radians(th)
    return X0 + s*np.cos(t) - n*np.sin(t), Y0 + s*np.sin(t) + n*np.cos(t)

def turn_back(u, v, rot, th=30.0):
    """velocity components of a turned run in the frame of the logical grid"""
    if rot == "rot90":
        return v, -u
    if rot == "rot":
        t = np.radians(th)
        return u*np.cos(t) + v*np.sin(t), -u*np.sin(t) + v*np.cos(t)
    return u, v

def box_phys(box, rot, X0=0.0, Y0=0.0, th=30.0):
    """F 72 box (xs xe ys ye) of the logical frame in the physical frame of a turned grid (rot90 exact; rot: the
    enclosing box is not the same set of cells, so the turned tests use initial states without boxes)"""
    xs, xe, ys, ye = box
    if rot == "rot90":
        return (X0-ye, X0-ys, Y0+xs, Y0+xe)
    return box

def run5(d, control, ctrl, grid=None, np_=1, steps=None, env_extra=None, files=None):
    """as p1lib.run, with the AMR dump of nhflow_amr / sflow_amr in d/amr (REEF3D_AMR_DUMP) and extra files"""
    if os.path.exists(d):
        shutil.rmtree(d)
    os.makedirs(os.path.join(d, "amr"))
    open(os.path.join(d, "control.txt"), "w").write(control)
    if steps is not None:
        ctrl = re.sub(r"^N 4[56] .*$", "", ctrl, flags=re.M) + "\nN 45 %d\n" % steps
    open(os.path.join(d, "ctrl.txt"), "w").write(ctrl)
    if grid is not None:
        write_grid(os.path.join(d, "curvgrid.dat"), *grid)
    for name, text in (files or {}).items():
        open(os.path.join(d, name), "w").write(text)
    with open(os.path.join(d, "dive.log"), "w") as log:
        r = subprocess.run([DIVEMESH], cwd=d, stdout=log, stderr=subprocess.STDOUT)
    if r.returncode != 0:
        return -100
    env = dict(os.environ)
    env["REEF3D_REGRESSION_DIR"] = os.path.join(d, "dump")
    env["REEF3D_AMR_DUMP"] = os.path.join(d, "amr")
    if env_extra:
        env.update(env_extra)
    cmd = [REEF3D] if np_ == 1 else ["mpirun", "--allow-run-as-root", "--oversubscribe", "-np", str(np_), REEF3D]
    with open(os.path.join(d, "reef.log"), "w") as log:
        r = subprocess.run(cmd, cwd=d, stdout=log, stderr=subprocess.STDOUT, env=env)
    return r.returncode

def amr_grids2(d, np_=1):
    """the grids of the AMR dump of all ranks (NHFLOW or SFLOW): list of dicts box (lev, blk, I0, J0, nx, ny), the 2D
    fields eta, WL, bed (nx, ny), and U, V (SFLOW: 2D; NHFLOW: the bottom layer)"""
    out = []
    for r in range(np_):
        st = read_state(os.path.join(d, "amr", "amr_r%d.bin" % r))
        ids = sorted({int(k[1:].split("_")[0]) for k in st if k.startswith("g")})
        for g in ids:
            lv, bk, I0, J0, nx, ny, nz = [int(v) for v in st["g%d_box" % g]]
            G = {"rank": r, "id": g, "lev": lv, "blk": bk, "I0": I0, "J0": J0, "nx": nx, "ny": ny, "count": st["count"], "simtime": st["simtime"]}
            for f in ("eta", "WL", "bed"):
                G[f] = st["g%d_%s" % (g, f)].reshape(nx, ny)
            for f in ("U", "V"):
                G[f] = st["g%d_%s" % (g, f)].reshape(nx, ny, nz)[:, :, 0]
            out.append(G)
    return out
