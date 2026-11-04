# the loaders of the Weber junction runs (weber_case.py: the AMR dump of the blocks; weber_cart.py: the VTP of the
# Cartesian reference), in the coordinates of Huang et al. (2002): x* = (Lu - x)/W positive upstream, the branch at
# x* -1 .. 0, y* = (W - y)/W from the branch-side wall
import glob, json, os, re, sys, numpy as np
T = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(T, "../sediment_curvilinear"))
from p5lib import amr_grids2, read_state
FT = 0.3048
W, Lm, Lb, Lu = 3*FT, 72*FT, 12*FT, 18*FT
U0 = 0.628      # the mean velocity of the main channel downstream of the junction (Weber et al. 2001: Q 0.170 m3/s, 0.296 m)

def read_blocks(path):
    """blocks.dat: the node arrays (ni + 1, nj + 1) of x and y per block"""
    t = open(path).read().split()
    k, out = 1, []
    for _ in range(int(t[0])):
        ni, nj = int(t[k]), int(t[k+1]); k += 2
        a = np.array(t[k:k+3*(ni+1)*(nj+1)], float).reshape(nj+1, ni+1, 3); k += 3*(ni+1)*(nj+1)
        out.append((a[:, :, 0].T, a[:, :, 1].T))
    return out

def load(d):
    """cell centres (x*, y* of Huang et al.), eta, depth and the velocity along x*, y* of all level-0 cells (the AMR
    dump has the Cartesian velocities); the cells from blocks.dat, block 0's from the native grid less its ring"""
    nodes = read_blocks(os.path.join(d, "blocks.dat"))
    ring = [0, 0, 0, 0]
    for l in open(os.path.join(d, "reef.log"), errors="replace"):
        m = re.search(r"ring (\d+) (\d+) (\d+) (\d+)", l)
        if m:
            ring = [int(v) for v in m.groups()]
    npr = len(glob.glob(os.path.join(d, "amr", "amr_r*.bin")))
    X, Y, Hd, U, V = [], [], [], [], []
    for g in grids_nz(d, npr):
        if g["lev"] != 0:
            continue
        # block 0: its cells from (I0, J0) of the native grid; the blocks 1..: their own cells (I0, J0 place the level-0
        # grid of a block in the index space of the hierarchy)
        o = (g["I0"] - ring[0], g["J0"] - ring[2]) if g["blk"] == 0 else (0, 0)
        nx_, ny_ = nodes[g["blk"]][0].shape[0]-1, nodes[g["blk"]][0].shape[1]-1
        I = o[0] + np.arange(g["nx"])[:, None] + 0*np.arange(g["ny"])[None, :]
        J = o[1] + 0*np.arange(g["nx"])[:, None] + np.arange(g["ny"])[None, :]
        keep = (I >= 0) & (I < nx_) & (J >= 0) & (J < ny_)
        x, y = nodes[g["blk"]]
        cx = 0.25*(x[:-1, :-1] + x[1:, :-1] + x[:-1, 1:] + x[1:, 1:])
        cy = 0.25*(y[:-1, :-1] + y[1:, :-1] + y[:-1, 1:] + y[1:, 1:])
        X.append(cx[I[keep], J[keep]]); Y.append(cy[I[keep], J[keep]])
        Hd.append(g["WL"][keep]); U.append(g["U"][keep]); V.append(g["V"][keep])
    X, Y, Hd, U, V = (np.concatenate(q) for q in (X, Y, Hd, U, V))
    dx = np.hypot(*(np.diff(nodes[0][0][:2, 0]), np.diff(nodes[0][1][:2, 0])))[0]
    xs, ys = (Lu - X)/W, (W - Y)/W                 # Huang's x*, y*
    # u*, v* along x*, y* (x* points upstream, y* away from the branch); NHFLOW: per sigma layer (u3, v3, the layers
    # of equal thickness from the bed) and their mean
    return dict(x=xs, y=ys, h=Hd, u=-U.mean(1), v=-V.mean(1), u3=-U, v3=-V, nz=U.shape[1], dx=dx/W, simtime=g["simtime"])

def grids_nz(d, npr):
    """the grids of the AMR dump with all layers: U, V (nx, ny, nz); SFLOW nz 1"""
    out = []
    for r in range(npr):
        st = read_state(os.path.join(d, "amr", "amr_r%d.bin" % r))
        ids = sorted({int(k[1:].split("_")[0]) for k in st if k.startswith("g")})
        for g in ids:
            lv, bk, I0, J0, nx, ny, nz = [int(v) for v in st["g%d_box" % g]]
            G = {"lev": lv, "blk": bk, "I0": I0, "J0": J0, "nx": nx, "ny": ny, "simtime": st["simtime"]}
            for f in ("eta", "WL", "bed"):
                G[f] = st["g%d_%s" % (g, f)].reshape(nx, ny)
            for f in ("U", "V"):
                G[f] = st["g%d_%s" % (g, f)].reshape(nx, ny, nz)
            out.append(G)
    return out

def read_vtp(fn):
    """the points and the point data of a VTK PolyData file of REEF3D (raw appended data, UInt32 headers)"""
    buf = open(fn, "rb").read()
    k = buf.index(b"<AppendedData")
    head = buf[:k].decode(errors="replace")
    start = buf.index(b"_", k) + 1
    out = {}
    for m in re.finditer(r'<DataArray type="(\w+)"(?: Name="(\w+)")?(?: NumberOfComponents="(\d)")? format="appended" offset="(\d+)"', head):
        typ, name, nc, off = m.group(1), m.group(2) or "Points", int(m.group(3) or 1), int(m.group(4))
        n = int(np.frombuffer(buf, "<u4", 1, start+off)[0])
        dt = {"Float32": "<f4", "Float64": "<f8", "Int32": "<i4"}[typ]
        out[name] = np.frombuffer(buf, dt, n//np.dtype(dt).itemsize, start+off+4).reshape(-1, nc) if nc > 1 else \
            np.frombuffer(buf, dt, n//np.dtype(dt).itemsize, start+off+4)
    return out

def load_cart(d):
    """the Cartesian reference of weber_cart.py: the last VTP of the free surface (point data), the points inside
    the channels"""
    meta = json.load(open(os.path.join(d, "case.json")))
    last = sorted(glob.glob(os.path.join(d, "REEF3D_SFLOW_VTP_FSF", "*.pvtp")))[-1]
    pieces = re.findall(r'Source="([^"]+)"', open(last).read())
    X, Y, Hd, U, V, t = [], [], [], [], [], 0.0
    for pc in pieces:
        v = read_vtp(os.path.join(os.path.dirname(last), pc))
        X.append(v["Points"][:, 0]); Y.append(v["Points"][:, 1]); Hd.append(v["waterlevel"])
        U.append(v["velocity"][:, 0]); V.append(v["velocity"][:, 1])
    X, Y, Hd, U, V = (np.concatenate(q).astype(float) for q in (X, Y, Hd, U, V))
    t = float(re.search(rb'Name="TimeValue"[^>]*>\s*([0-9.eE+-]+)', open(os.path.join(os.path.dirname(last), pieces[0]), "rb").read()).group(1))
    e = 1e-6
    inside = (Y <= W+e) | ((X >= Lu-e) & (X <= Lu+W+e))
    X, Y, Hd, U, V = (q[inside] for q in (X, Y, Hd, U, V))
    xs, ys = (Lu - X)/W, (W - Y)/W
    return dict(x=xs, y=ys, h=Hd, u=-U, v=-V, dx=meta["dx"]/W, simtime=t, note=meta["note"])

def at(r, xs, ys):
    """nearest cell"""
    k = np.argmin((r["x"]-xs)**2 + (r["y"]-ys)**2)
    return k

def load_any(d):
    return load_cart(d) if os.path.exists(os.path.join(d, "case.json")) else load(d)
