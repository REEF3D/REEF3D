# Reader of a REEF3D LAGOON store (P 18, Zarr v3, format 0.5), numpy only:
#   s = lagoon("REEF3D_CFD.lagoon")
#   s.times                      output times
#   s.coords()                   x, y, z of the grid points of the CFD volume
#   s.volume("pressure", n)      field of output n on the whole grid, shape (z, y, x[, components])
#   s.body(n)                    vertices (moved by the body motion of output n) and triangles
#   s.amr_times()                output times of the AMR set (cfd_amr: the CFD mesh refinement)
#   s.amr(n, ["pressure", ...])  the grids of AMR output n (level 0 of every rank, then the
#                                patches): dicts with level, rank, x, y, z (nodes), the fields (z, y, x)
# Arrays: sharded (sharding_indexed: inner chunks, index at the end with crc32c), each inner chunk
# bytes little endian -> numcodecs.shuffle -> gzip, as REEF3D writes them (lagoon_store.cpp).
import json, os, zlib
import numpy as np

EMPTY = np.uint64(0xFFFFFFFFFFFFFFFF)


class zarr_array:
    def __init__(self, path):
        self.path = path
        m = json.load(open(os.path.join(path, "zarr.json")))
        self.shape = m["shape"]
        self.dtype = np.dtype(m["data_type"]).newbyteorder("<")
        self.chunk = m["chunk_grid"]["configuration"]["chunk_shape"]
        self.fill = m.get("fill_value")
        self.sharded = m["codecs"][0]["name"] == "sharding_indexed"
        self.inner_codecs = [x["name"] for x in m["codecs"]]
        if self.sharded:
            c = m["codecs"][0]["configuration"]
            self.inner = c["chunk_shape"]
            self.inner_codecs = [x["name"] for x in c["codecs"]]
            self.nin = [s//i for s, i in zip(self.chunk, self.inner)]

    def _key(self, cidx):
        return os.path.join(self.path, "c", *[str(i) for i in cidx])

    def _decode(self, buf, n):
        if "gzip" in self.inner_codecs:
            buf = zlib.decompress(buf, 31)
        if "numcodecs.shuffle" in self.inner_codecs:
            es = self.dtype.itemsize
            buf = np.frombuffer(buf, np.uint8).reshape(es, -1).T.tobytes()
        a = np.frombuffer(buf, self.dtype)
        return a[:n] if len(a) >= n else np.concatenate([a, np.full(n - len(a), np.nan if self.dtype.kind == "f" else 0, self.dtype)])

    def read_all(self):
        # unsharded small arrays (coordinates, time)
        out = np.full(self.shape, np.nan if self.dtype.kind == "f" else 0, self.dtype)
        nch = [-(-s//c) for s, c in zip(self.shape, self.chunk)]
        for idx in np.ndindex(*nch):
            f = self._key(idx)
            if not os.path.exists(f):
                continue
            # a chunk holds the full chunk shape (fill beyond the array), decoded by the codecs
            a = self._decode(open(f, "rb").read(), int(np.prod(self.chunk))).reshape(self.chunk)
            sl = tuple(slice(i*c, min((i+1)*c, s)) for i, c, s in zip(idx, self.chunk, self.shape))
            out[sl] = a[tuple(slice(0, x.stop - x.start) for x in sl)]
        return out

    def read_range(self, lo, hi):
        # rows lo to hi of the first dimension (sharded arrays, inner chunks over the first dimension only)
        rest = tuple(self.shape[1:])
        out = np.full((hi - lo,) + rest, np.nan if self.dtype.kind == "f" else 0, self.dtype)
        if hi <= lo:
            return out
        S, I = self.chunk[0], self.inner[0]
        ntot = int(np.prod(self.nin))
        n_inner = int(np.prod(self.inner))
        for sh in range(lo//S, (hi - 1)//S + 1):
            f = self._key([sh] + [0]*len(rest))
            if not os.path.exists(f):
                continue
            b = open(f, "rb").read()
            idx = np.frombuffer(b[-(16*ntot+4):-4], "<u8").reshape(ntot, 2)
            for k, mi in enumerate(np.ndindex(*self.nin)):
                r0 = sh*S + mi[0]*I
                a0, a1 = max(lo, r0), min(hi, r0 + I, self.shape[0])
                if a1 <= a0 or any(mi[1:]):
                    continue
                off, nb = idx[k]
                if off == EMPTY:
                    continue
                a = self._decode(b[int(off):int(off)+int(nb)], n_inner).reshape(self.inner)
                out[a0-lo:a1-lo] = a[a0-r0:a1-r0][(slice(None),) + tuple(slice(0, r) for r in rest)]
        return out

    def read_time(self, t):
        # all of output t (first dimension); sharded arrays
        ct = t//self.chunk[0]
        f = self._key([ct] + [0]*(len(self.shape)-1))
        b = open(f, "rb").read()
        ntot = int(np.prod(self.nin))
        idx = np.frombuffer(b[-(16*ntot+4):-4], "<u8").reshape(ntot, 2)
        out = np.full(self.shape[1:], np.nan, self.dtype)
        tin = t - ct*self.chunk[0]
        n_inner = int(np.prod(self.inner))
        for k, mi in enumerate(np.ndindex(*self.nin)):
            if mi[0] != tin // self.inner[0]:
                continue
            off, nb = idx[k]
            if off == EMPTY:
                continue
            a = self._decode(b[int(off):int(off)+int(nb)], n_inner).reshape(self.inner)
            a = a[tin % self.inner[0]]
            # inner chunks at the end of a dimension reach past the array: clipped
            sl = tuple(slice(i*s, min((i+1)*s, n)) for i, s, n in zip(mi[1:], self.inner[1:], self.shape[1:]))
            out[sl] = a[tuple(slice(0, x.stop - x.start) for x in sl)]
        return out


class lagoon:
    def __init__(self, store):
        self.store = store
        v = os.path.join(store, "volume")
        self.times = zarr_array(os.path.join(v, "time")).read_all()
        self.meta = json.load(open(os.path.join(v, "zarr.json")))["attributes"]["lagoon"]
        self.blocks = []
        bdir = os.path.join(v, "blocks")
        for b in sorted(d for d in os.listdir(bdir) if d.startswith("b")):
            a = json.load(open(os.path.join(bdir, b, "zarr.json")))["attributes"]["lagoon"]
            self.blocks.append((os.path.join(bdir, b), a["start"], a["size"]))

    def coords(self):
        v = os.path.join(self.store, "volume")
        return tuple(zarr_array(os.path.join(v, c)).read_all() for c in ("x", "y", "z"))

    def volume(self, var, n):
        sh = self.meta["shape"]
        comp = self.meta["variables"][var].get("components", 1)
        full = np.full((sh["level"], sh["y"], sh["x"]) + ((comp,) if comp > 1 else ()), np.nan, np.float32)
        for path, start, size in self.blocks:
            a = zarr_array(os.path.join(path, var)).read_time(n)
            full[start[2]:start[2]+size[2], start[1]:start[1]+size[1], start[0]:start[0]+size[0]] = a
        return full

    def body(self, n, ib=0):
        d = os.path.join(self.store, "bodies", "cfd_body", f"body_{ib}")
        V = zarr_array(os.path.join(d, "vertices")).read_all()
        T = zarr_array(os.path.join(d, "triangles")).read_all()
        tr = zarr_array(os.path.join(d, "translation")).read_all()
        q = zarr_array(os.path.join(d, "rotation")).read_all()
        t = zarr_array(os.path.join(self.store, "bodies", "cfd_body", "time")).read_all()
        # the body at output n: rotation (quaternion w, x, y, z) and translation of the nearest body time
        k = int(np.argmin(np.abs(t - self.times[n])))
        w, qx, qy, qz = q[k]
        Rm = np.array([[1-2*(qy*qy+qz*qz), 2*(qx*qy-qz*w), 2*(qx*qz+qy*w)],
                       [2*(qx*qy+qz*w), 1-2*(qx*qx+qz*qz), 2*(qy*qz-qx*w)],
                       [2*(qx*qz-qy*w), 2*(qy*qz+qx*w), 1-2*(qx*qx+qy*qy)]])
        return V @ Rm.T + tr[k], T

    # AMR set (P 18 with mesh refinement): flat arrays, every grid of every output one after another
    def _amr(self, name, key="cfd_amr"):
        return zarr_array(os.path.join(self.store, "amr", key, name))

    def amr_times(self, key="cfd_amr"):
        return self._amr("time", key).read_all()

    def amr(self, n, fields, key="cfd_amr"):
        ge = self._amr("grid_end", key).read_all()
        g0, g1 = (int(ge[n-1]) if n > 0 else 0), int(ge[n])
        lo = max(g0 - 1, 0)
        lev = self._amr("grids/level", key).read_range(g0, g1)
        rank = self._amr("grids/rank", key).read_range(g0, g1)
        ends = {e: self._amr(f"grids/{e}_end", key).read_range(lo, g1).astype(np.int64)
                for e in ("x", "y", "z", "cell")}
        if g0 == 0:
            ends = {e: np.concatenate([[0], v]) for e, v in ends.items()}
        # node and cell ranges of the grids of this output
        nodes = {d: self._amr(d, key).read_range(int(ends[d][0]), int(ends[d][-1])) for d in ("x", "y", "z")}
        c0 = int(ends["cell"][0])
        vals = {f: self._amr(f"cell_data/{f}", key).read_range(c0, int(ends["cell"][-1])) for f in fields}
        out = []
        for g in range(g1 - g0):
            gr = {"level": int(lev[g]), "rank": int(rank[g])}
            for d in ("x", "y", "z"):
                gr[d] = nodes[d][ends[d][g]-ends[d][0]:ends[d][g+1]-ends[d][0]]
            sh = (len(gr["z"]) - 1, len(gr["y"]) - 1, len(gr["x"]) - 1)
            for f in fields:
                gr[f] = vals[f][ends["cell"][g]-c0:ends["cell"][g+1]-c0].reshape(sh)
            out.append(gr)
        return out
