# Minimal reader of the REEF3D::CFD VTU output (appended raw binary, one .vtu per rank and a .pvtu
# per output step), numpy only:
#   steps(dir)          -> list of (step number, time) from the .pvtu files
#   read_step(dir, n)   -> dict of the point arrays of output n over all ranks, plus "points"
import glob, os, re
import numpy as np


def read_piece(fn):
    b = open(fn, "rb").read()
    k = b.index(b"<AppendedData")
    head = b[:k].decode("latin1")
    base = b.index(b"_", k) + 1
    out = {}
    for m in re.finditer(r'<DataArray type="Float32"( Name="([^"]+)")?( NumberOfComponents="(\d)")? '
                         r'format="appended" offset="(\d+)"', head):
        name = m.group(2) or "points"
        nc = int(m.group(4) or 1)
        off = int(m.group(5))
        n = np.frombuffer(b, dtype="<i4", count=1, offset=base+off)[0]
        a = np.frombuffer(b, dtype="<f4", count=n//4, offset=base+off+4)
        out[name] = a.reshape(-1, nc) if nc > 1 else a
    return out


def steps(d="REEF3D_CFD_VTU"):
    out = []
    for f in sorted(glob.glob(os.path.join(d, "REEF3D-CFD-*.pvtu"))):
        n = int(re.search(r"REEF3D-CFD-(\d+)\.pvtu", f).group(1))
        m = re.search(r'Name="TimeValue"[^>]*>\s*([-+0-9.eE]+)', open(f).read())
        out.append((n, float(m.group(1)) if m else float(n)))
    return out


def read_step(d, n):
    parts = [read_piece(f) for f in sorted(glob.glob(os.path.join(d, f"REEF3D-CFD-{n:08d}-*.vtu")))]
    return {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
