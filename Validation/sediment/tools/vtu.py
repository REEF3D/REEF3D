import numpy as np, re
def read_vtu(fn):
    raw=open(fn,'rb').read(); i=raw.index(b'<AppendedData'); j=raw.index(b'_',i)+1
    head=raw[:i].decode('latin1')
    out={}; 
    for m in re.finditer(r'<DataArray([^>]*)/?>',head):
        a=m.group(1); off=re.search(r'offset="(\d+)"',a)
        if off is None: continue
        name=re.search(r'Name="([^"]*)"',a); typ=re.search(r'type="([^"]*)"',a).group(1)
        nc=re.search(r'NumberOfComponents="(\d+)"',a); nc=int(nc.group(1)) if nc else 1
        off=int(off.group(1)); n=np.frombuffer(raw[j+off:j+off+4],dtype=np.int32)[0]
        dt={'Float32':np.float32,'Float64':np.float64,'Int32':np.int32,'Int64':np.int64,'UInt8':np.uint8}[typ]
        data=np.frombuffer(raw[j+off+4:j+off+4+n],dtype=dt)
        key=name.group(1) if name else 'points'
        if key in out: key=key+'_2'
        out[key]=data.reshape(-1,nc) if nc>1 else data
    return out
