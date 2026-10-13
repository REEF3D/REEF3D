# helpers for the multi-block NHFLOW tests (G 50 1 with DIVEMesh R 20 1, phase P4): runs go to ./runs
# binaries: environment REEF3D_BIN, DIVEMESH_BIN (DIVEMesh with R 20, the multi-block grid)
import os, sys, re, glob, shutil, subprocess, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../sflow_multiblock"))
from p2lib import *
from p2lib import write_blocks
S4 = os.path.dirname(os.path.abspath(__file__))

def run4(d, control, ctrl, blocks, np_=1, steps=None, env_extra=None):
    """DIVEMesh and REEF3D in a fresh folder d; the grids of nhflow_amr (REEF3D_AMR_DUMP) go to d/amr"""
    d = os.path.abspath(d)
    if os.path.exists(d):
        shutil.rmtree(d)
    os.makedirs(os.path.join(d, "amr"))
    open(os.path.join(d, "control.txt"), "w").write(control)
    if steps is not None:
        ctrl = re.sub(r"^N 4[56] .*$", "", ctrl, flags=re.M) + "\nN 45 %d\n" % steps
    open(os.path.join(d, "ctrl.txt"), "w").write(ctrl)
    write_blocks(os.path.join(d, "blocks.dat"), blocks)
    with open(os.path.join(d, "dive.log"), "w") as log:
        r = subprocess.run([DM], cwd=d, stdout=log, stderr=subprocess.STDOUT)
    if r.returncode != 0:
        return -100
    env = dict(os.environ)
    env["REEF3D_REGRESSION_DIR"] = os.path.join(d, "dump")
    env["REEF3D_AMR_DUMP"] = os.path.join(d, "amr")
    if env_extra:
        env.update(env_extra)
    mpi = os.environ.get("REEF3D_MPIRUN", "mpirun --allow-run-as-root --oversubscribe").split()
    cmd = [REEF3D] if np_ == 1 else mpi + ["-np", str(np_), REEF3D]
    with open(os.path.join(d, "reef.log"), "w") as log:
        r = subprocess.run(cmd, cwd=d, stdout=log, stderr=subprocess.STDOUT, env=env)
    return r.returncode

def amr_grids(d, np_=1):
    """the grids of all ranks: list of dicts box (lev, blk, I0, J0, nx, ny, nz), eta, WL, bed (nx, ny),
    U, V, W, UH, VH, WH (nx, ny, nz), P (nx, ny, nz+1)"""
    out = []
    for r in range(np_):
        st = read_state(os.path.join(d, "amr", "amr_r%d.bin" % r))
        ids = sorted({int(k[1:].split("_")[0]) for k in st if k.startswith("g")})
        for g in ids:
            b = st["g%d_box" % g]
            lv, bk, I0, J0, nx, ny, nz = [int(v) for v in b]
            G = {"rank": r, "id": g, "lev": lv, "blk": bk, "I0": I0, "J0": J0, "nx": nx, "ny": ny, "nz": nz, "count": st["count"], "simtime": st["simtime"]}
            for f in ("eta", "WL", "bed"):
                G[f] = st["g%d_%s" % (g, f)].reshape(nx, ny)
            for f in ("U", "V", "W", "UH", "VH", "WH"):
                G[f] = st["g%d_%s" % (g, f)].reshape(nx, ny, nz)
            G["P"] = st["g%d_P" % g].reshape(nx, ny, nz+1)
            out.append(G)
    return out

def level0(grids, gnx, gny, ring=(0, 0, 0, 0)):
    """the native level-0 grid assembled from the ranks (block 0 with its ring), fields (GNX, GNY, ...)"""
    out = {}
    for G in grids:
        if G["id"] != 0:
            continue
        for f in ("eta", "WL", "bed", "U", "V", "W", "UH", "VH", "WH", "P"):
            a = G[f]
            if f not in out:
                out[f] = np.zeros((gnx, gny) + a.shape[2:])
            out[f][G["I0"]:G["I0"]+G["nx"], G["J0"]:G["J0"]+G["ny"]] = a
    return out

FIELDS = ("eta", "WL", "U", "V", "W", "P")
