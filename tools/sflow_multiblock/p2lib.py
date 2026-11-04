# helpers for the multi-block SFLOW tests (G 50 1 with DIVEMesh R 20 1): runs go to ./runs
# binaries: environment REEF3D_BIN, DIVEMESH_BIN (DIVEMesh with R 20, the multi-block grid)
import os, sys, re, shutil, subprocess, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../sflow_curvilinear"))
from p1lib import *
S2 = os.path.dirname(os.path.abspath(__file__))
DM = DIVEMESH

def write_blocks(path, blocks):
    """blocks.dat: number of blocks, per block "ni nj" and the nodes "x y zb" (i fastest);
    blocks: list of (X,Y,Z) arrays (ni+1, nj+1)"""
    with open(path, "w") as f:
        f.write("%d\n" % len(blocks))
        for X, Y, Z in blocks:
            ni, nj = X.shape[0]-1, X.shape[1]-1
            f.write("%d %d\n" % (ni, nj))
            for j in range(nj+1):
                for i in range(ni+1):
                    f.write("%.17g %.17g %.17g\n" % (X[i, j], Y[i, j], Z[i, j]))

def run_mb(d, control, ctrl, blocks, np_=1, steps=None, reef=None):
    """DIVEMesh and REEF3D in a fresh folder d; the regression dump of REEF3D goes to d/dump"""
    if os.path.exists(d):
        shutil.rmtree(d)
    os.makedirs(d)
    open(os.path.join(d, "control.txt"), "w").write(control)
    if steps is not None:
        ctrl = re.sub(r"^N 4[56] .*$", "", ctrl, flags=re.M) + "\nN 45 %d\n" % steps
    open(os.path.join(d, "ctrl.txt"), "w").write(ctrl)
    write_blocks(os.path.join(d, "blocks.dat"), blocks)
    with open(os.path.join(d, "dive.log"), "w") as log:
        r = subprocess.run([DM], cwd=d, stdout=log, stderr=subprocess.STDOUT)
    if r.returncode != 0:
        return -100
    env = dict(os.environ); env["REEF3D_REGRESSION_DIR"] = os.path.join(d, "dump")
    mpi = os.environ.get("REEF3D_MPIRUN", "mpirun --allow-run-as-root --oversubscribe").split()
    cmd = [reef or REEF3D] if np_ == 1 else mpi + ["-np", str(np_), reef or REEF3D]
    with open(os.path.join(d, "reef.log"), "w") as log:
        r = subprocess.run(cmd, cwd=d, stdout=log, stderr=subprocess.STDOUT, env=env)
    return r.returncode
