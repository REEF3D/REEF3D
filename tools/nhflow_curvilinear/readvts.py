import re, numpy as np, sys
def readvts(fn):
    s=open(fn).read()
    ext=[int(v) for v in re.search(r'WholeExtent="([^"]+)"',s).group(1).split()]
    nx,ny=ext[1],ext[3]
    out={}
    for m in re.finditer(r'<DataArray type="Float64" Name="(\w+)"[^>]*format="ascii">(.*?)</DataArray>',s,re.S):
        a=np.array(m.group(2).split(),float)
        if a.size==nx*ny: out[m.group(1)]=a.reshape(ny,nx).T
    return out,nx,ny
