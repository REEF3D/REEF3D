# AMR VTS output (ascii) of REEF3D::SFLOW: cell arrays, cell centres and level
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

def readvts_full(fn):
    """cell arrays, cell centres (mean of the 4 nodes) and the level of an AMR VTS file"""
    s=open(fn).read()
    out,nx,ny=readvts(fn)
    m=re.search(r'<Points>\s*<DataArray[^>]*>(.*?)</DataArray>',s,re.S)
    P=np.array(m.group(1).split(),float).reshape(ny+1,nx+1,3)
    xc=0.25*(P[:-1,:-1,0]+P[1:,:-1,0]+P[:-1,1:,0]+P[1:,1:,0]).T
    yc=0.25*(P[:-1,:-1,1]+P[1:,:-1,1]+P[:-1,1:,1]+P[1:,1:,1]).T
    lv=re.search(r'Name="level"[^>]*>\s*(-?\d+)',s)
    return out,xc,yc,(int(lv.group(1)) if lv else -1)
