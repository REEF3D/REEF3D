# Architect: Hans Bihs
"""Readers for REEF3D output used by the turbulence test cases.

read_vtk(path)        VTK XML file (vtu/vtp) with raw appended data, as written by REEF3D
                      -> dict: 'points' (n,3) and every DataArray by name (point or cell data)
read_state(path)      regression_dump state_<count>_r<rank>.bin -> dict name -> array
"""
import re, struct
import numpy as np

_TYPES = {"Float32": np.float32, "Float64": np.float64, "Int32": np.int32, "Int64": np.int64,
          "UInt8": np.uint8, "Int8": np.int8, "UInt32": np.uint32, "UInt64": np.uint64}


def read_vtk(path):
    raw = open(path, "rb").read()
    m = re.search(rb"<AppendedData[^>]*>\s*_", raw)
    if not m:
        raise ValueError("no raw appended data in " + path)
    base = m.end()
    head = raw[:m.start()].decode("latin-1")
    hdr = "UInt64" if 'header_type="UInt64"' in head else "UInt32"
    hfmt, hsize = ("<Q", 8) if hdr == "UInt64" else ("<I", 4)
    out = {}
    # points: first DataArray inside <Points>
    pm = re.search(r"<Points>\s*<DataArray([^>]*)>", head)
    arrays = re.findall(r"<DataArray([^>]*)/?>", head)
    def parse(attrs):
        d = dict(re.findall(r'(\w+)="([^"]*)"', attrs))
        return d
    for a in arrays:
        d = parse(a)
        if d.get("format") != "appended":
            continue
        off = int(d["offset"])
        nbytes, = struct.unpack_from(hfmt, raw, base + off)
        dt = _TYPES[d["type"]]
        arr = np.frombuffer(raw, dtype=dt, count=nbytes // np.dtype(dt).itemsize, offset=base + off + hsize)
        nc = int(d.get("NumberOfComponents", "1"))
        if nc > 1:
            arr = arr.reshape(-1, nc)
        name = d.get("Name")
        if pm and a == pm.group(1):
            name = "points"
        if name is None:
            continue
        out[name] = arr.astype(np.float64) if arr.dtype.kind == "f" else arr
    return out


def read_state(path):
    buf = open(path, "rb").read()
    if buf[:8] != b"R3DREG01":
        raise ValueError("bad regression state file " + path)
    pos = 8
    rank, size, count = struct.unpack_from("<iii", buf, pos); pos += 12
    simtime, = struct.unpack_from("<d", buf, pos); pos += 8
    nf, = struct.unpack_from("<i", buf, pos); pos += 4
    out = {"_count": count, "_simtime": simtime}
    for _ in range(nf):
        name = buf[pos:pos + 16].split(b"\0")[0].decode(); pos += 16
        n, = struct.unpack_from("<q", buf, pos); pos += 8
        out[name] = np.frombuffer(buf, dtype="<f8", count=n, offset=pos).copy(); pos += 8 * n
    return out
