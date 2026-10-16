import struct, numpy as np, sys
def read_curv(fn):
    b = open(fn, "rb").read()
    assert b[:8] == b"DMGEOM02", b[:8]
    pos = 16
    while pos + 12 <= len(b):
        tag = b[pos:pos+4].decode(); n, = struct.unpack_from("<q", b, pos+4); pos += 12
        if tag == "CURV":
            q = pos
            nint, = struct.unpack_from("<i", b, q); q += 4
            iv = struct.unpack_from("<%di" % nint, b, q); q += 4*nint
            nd, = struct.unpack_from("<i", b, q); q += 4
            par = struct.unpack_from("<%dd" % nd, b, q); q += 8*nd
            ni, nj = iv[1], iv[2]; nn = (ni+1)*(nj+1)
            arr = np.frombuffer(b[q:q+24*nn], dtype="<f8").reshape(3, nj+1, ni+1)
            return ni, nj, par, arr[0].T.copy(), arr[1].T.copy(), arr[2].T.copy()
        if tag == "END ":
            break
        pos += n
    raise ValueError("no CURV")
if __name__ == "__main__":
    ni, nj, par, X, Y, Z = read_curv(sys.argv[1])
    print(ni, nj, par, X.shape, X.min(), X.max(), Y.min(), Y.max(), Z.min(), Z.max())
