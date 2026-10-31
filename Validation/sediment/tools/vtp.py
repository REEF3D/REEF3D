import numpy as np, re, glob, os
def read_vtp(fn):
    """REEF3D appended-raw VTP (serial piece): returns points (N,3) and dict of point data arrays"""
    raw=open(fn,'rb').read()
    i=raw.index(b'<AppendedData')
    j=raw.index(b'_',i)+1
    head=raw[:i].decode('latin1')
    npts=int(re.search(r'NumberOfPoints="(\d+)"',head).group(1))
    arrays=[]
    for m in re.finditer(r'<DataArray([^>]*)>',head):
        a=m.group(1)
        name=re.search(r'Name="([^"]*)"',a); off=re.search(r'offset="(\d+)"',a)
        typ=re.search(r'type="([^"]*)"',a).group(1)
        nc=re.search(r'NumberOfComponents="(\d+)"',a)
        if off is None: continue
        arrays.append((name.group(1) if name else '?',int(off.group(1)),typ,int(nc.group(1)) if nc else 1))
    out={}
    for name,off,typ,nc in arrays:
        n=np.frombuffer(raw[j+off:j+off+4],dtype=np.int32)[0]
        dt={'Float32':np.float32,'Float64':np.float64,'Int32':np.int32,'Int64':np.int64,'UInt8':np.uint8}.get(typ,np.float32)
        data=np.frombuffer(raw[j+off+4:j+off+4+n],dtype=dt)
        out[name]=data.reshape(-1,nc) if nc>1 else data
    return out
def last_vtp(d):
    fs=sorted(glob.glob(os.path.join(d,'*.vtp')),key=lambda f:int(re.findall(r'(\d+)',os.path.basename(f))[-1]) if re.findall(r'(\d+)',os.path.basename(f)) else 0)
    return fs
