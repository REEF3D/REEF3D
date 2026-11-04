# helpers for the curvilinear SFLOW tests (G 50): runs go to ./runs
import os, subprocess, struct, glob, re, shutil, sys
import numpy as np

S = os.path.dirname(os.path.abspath(__file__))
# binaries: environment REEF3D_BIN, DIVEMESH_BIN (default: REEF3D of this repository, DiveMESH in PATH)
DIVEMESH = os.environ.get("DIVEMESH_BIN", shutil.which("DiveMESH") or "DiveMESH")
REEF3D = os.environ.get("REEF3D_BIN", os.path.join(S, "../../bin/REEF3D"))

def write_grid(path, X, Y, Z):
    """X,Y,Z: arrays (ni+1, nj+1) of node coordinates and bed levels"""
    ni, nj = X.shape[0]-1, X.shape[1]-1
    with open(path, "w") as f:
        f.write("%d %d\n" % (ni, nj))
        for j in range(nj+1):
            for i in range(ni+1):
                f.write("%.17g %.17g %.17g\n" % (X[i, j], Y[i, j], Z[i, j]))

def run(d, control, ctrl, grid=None, np_=1, steps=None, env_extra=None):
    if os.path.exists(d):
        shutil.rmtree(d)
    os.makedirs(d)
    open(os.path.join(d, "control.txt"), "w").write(control)
    if steps is not None:
        ctrl = re.sub(r"^N 4[56] .*$", "", ctrl, flags=re.M) + "\nN 45 %d\n" % steps
    open(os.path.join(d, "ctrl.txt"), "w").write(ctrl)
    if grid is not None:
        write_grid(os.path.join(d, "curvgrid.dat"), *grid)
    with open(os.path.join(d, "dive.log"), "w") as log:
        subprocess.run([DIVEMESH], cwd=d, stdout=log, stderr=subprocess.STDOUT, check=True)
    env = dict(os.environ)
    env["REEF3D_REGRESSION_DIR"] = os.path.join(d, "dump")
    if env_extra:
        env.update(env_extra)
    cmd = [REEF3D] if np_ == 1 else ["mpirun", "--allow-run-as-root", "--oversubscribe", "-np", str(np_), REEF3D]
    with open(os.path.join(d, "reef.log"), "w") as log:
        r = subprocess.run(cmd, cwd=d, stdout=log, stderr=subprocess.STDOUT, env=env)
    return r.returncode

def read_state(fn):
    buf = open(fn, "rb").read()
    assert buf[:8] == b"R3DREG01"
    pos = 8
    rank, size, count = struct.unpack_from("<iii", buf, pos); pos += 12
    simtime, = struct.unpack_from("<d", buf, pos); pos += 8
    nf, = struct.unpack_from("<i", buf, pos); pos += 4
    out = {"count": count, "simtime": simtime}
    for _ in range(nf):
        name = buf[pos:pos+16].split(b"\0")[0].decode(); pos += 16
        n, = struct.unpack_from("<q", buf, pos); pos += 8
        out[name] = np.frombuffer(buf[pos:pos+8*n], dtype="<f8").copy(); pos += 8*n
    return out

def final_state(d, rank=0):
    fs = glob.glob(os.path.join(d, "dump", "state_*_r%d.bin" % rank))
    c = max(int(re.search(r"state_(\d+)_r", f).group(1)) for f in fs)
    return read_state(os.path.join(d, "dump", "state_%08d_r%d.bin" % (c, rank)))

def steps(d, rank=0):
    rows = []
    for l in open(os.path.join(d, "dump", "steps_r%d.txt" % rank)):
        if l.startswith("#") or not l.strip():
            continue
        t = l.split()
        rows.append([float(t[0])] + [float.fromhex(x) for x in t[1:]])
    return np.array(rows)

def log_lines(d, pat):
    return [l.rstrip() for l in open(os.path.join(d, "reef.log"), errors="replace") if re.search(pat, l)]
