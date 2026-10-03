/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"fem_solid.h"
#include<cmath>
#include<cstdint>
#include<cstring>
#include<fstream>
#include<sstream>
#include<stdexcept>
#include<algorithm>
#include<array>

namespace
{
    typedef std::array<double,9> tri;   // x0 y0 z0 x1 y1 z1 x2 y2 z2

    std::vector<tri> read_stl(const std::string& file)
    {
        std::ifstream f(file.c_str(),std::ios::binary);
        if(!f)
        throw std::runtime_error("FEM: cannot open STL file "+file);

        f.seekg(0,std::ios::end);
        const std::streamoff size = f.tellg();
        f.seekg(0,std::ios::beg);

        std::vector<tri> T;

        // binary STL: 80 byte header, uint32 count, 50 bytes per facet
        if(size>=84)
        {
            char header[80];
            uint32_t n = 0;
            f.read(header,80);
            f.read(reinterpret_cast<char*>(&n),4);

            if(size==84+50*(std::streamoff)n)
            {
                T.resize(n);
                for(uint32_t t=0; t<n; ++t)
                {
                    float buf[12];
                    uint16_t attr;
                    f.read(reinterpret_cast<char*>(buf),48);
                    f.read(reinterpret_cast<char*>(&attr),2);
                    for(int q=0; q<9; ++q)
                    T[t][q] = buf[3+q];
                }
                return T;
            }
        }

        // ASCII STL
        f.clear();
        f.seekg(0,std::ios::beg);
        std::string word;
        tri cur;
        int nv = 0;
        while(f>>word)
        {
            if(word=="vertex")
            {
                f>>cur[3*nv]>>cur[3*nv+1]>>cur[3*nv+2];
                ++nv;
                if(nv==3)
                {
                    T.push_back(cur);
                    nv = 0;
                }
            }
        }

        if(T.empty())
        throw std::runtime_error("FEM: no triangles found in STL file "+file);

        return T;
    }

    void stl_bounds(const std::vector<tri>& T,double* lo,double* hi)
    {
        for(int d=0; d<3; ++d) {lo[d] = 1.0e300; hi[d] = -1.0e300;}
        for(const tri& t : T)
        for(int v=0; v<3; ++v)
        for(int d=0; d<3; ++d)
        {
            lo[d] = std::min(lo[d],t[3*v+d]);
            hi[d] = std::max(hi[d],t[3*v+d]);
        }
    }
}

int fem_solid::voxel(int ix,int iy,int iz) const
{
    if(ix<0 || iy<0 || iz<0 || ix>=nx || iy>=ny || iz>=nz)
    return -1;
    return ix + nx*(iy + ny*iz);
}

void fem_solid::voxelise()
{
    // lattice extent: union of all added shapes
    double lo[3] = {1.0e300,1.0e300,1.0e300};
    double hi[3] = {-1.0e300,-1.0e300,-1.0e300};
    bool any = false;

    std::vector<std::vector<tri>> stltri(stls.size());

    for(size_t s=0; s<stls.size(); ++s)
    stltri[s] = read_stl(stls[s].file);

    for(const shape_cmd& c : shapes)
    {
        double l[3], h[3];
        bool add = false;

        if(c.type==0 && !boxes[c.idx].remove)
        {
            const shape_box& b = boxes[c.idx];
            l[0]=b.x0; l[1]=b.y0; l[2]=b.z0; h[0]=b.x1; h[1]=b.y1; h[2]=b.z1;
            add = true;
        }
        if(c.type==1 && !stls[c.idx].remove)
        {
            stl_bounds(stltri[c.idx],l,h);
            add = true;
        }

        if(add)
        {
            any = true;
            for(int d=0; d<3; ++d) {lo[d] = std::min(lo[d],l[d]); hi[d] = std::max(hi[d],h[d]);}
        }
    }

    if(!any)
    throw std::runtime_error("FEM: no geometry (keywords 'box' or 'stl')");

    // origin: given lattice origin shifted by whole cells to just below the geometry
    const double hh[3] = {hx,hy,hz};
    double org[3] = {ox,oy,oz};
    int nn[3];
    for(int d=0; d<3; ++d)
    {
        const double shift = std::floor((lo[d]-org[d])/hh[d] + 1.0e-9);
        org[d] += shift*hh[d];
        nn[d] = std::max(1,(int)std::ceil((hi[d]-org[d])/hh[d] - 1.0e-9));
    }
    ox=org[0]; oy=org[1]; oz=org[2];
    nx=nn[0]; ny=nn[1]; nz=nn[2];

    if((double)nx*(double)ny*(double)nz > 2.0e8)
    throw std::runtime_error("FEM: voxel lattice too large");

    vox.assign((size_t)nx*ny*nz,-1);

    auto matindex = [&](int id)->int
    {
        for(size_t n=0; n<mats.size(); ++n)
        if(mats[n].id==id) return (int)n;
        throw std::runtime_error("FEM: material "+std::to_string(id)+" not defined");
    };

    for(const shape_cmd& c : shapes)
    {
        if(c.type==0)
        {
            const shape_box& b = boxes[c.idx];
            const int mi = b.remove ? -1 : matindex(b.mat);

            for(int k=0; k<nz; ++k)
            for(int j=0; j<ny; ++j)
            for(int i=0; i<nx; ++i)
            {
                const double xc = ox+(i+0.5)*hx, yc = oy+(j+0.5)*hy, zc = oz+(k+0.5)*hz;
                if(xc>b.x0 && xc<b.x1 && yc>b.y0 && yc<b.y1 && zc>b.z0 && zc<b.z1)
                vox[voxel(i,j,k)] = mi;
            }
        }

        if(c.type==1)
        {
            const shape_stl& s = stls[c.idx];
            const int mi = s.remove ? -1 : matindex(s.mat);
            const std::vector<tri>& T = stltri[c.idx];

            // ray parity along +z through every voxel column
            std::vector<double> zhit;
            for(int j=0; j<ny; ++j)
            for(int i=0; i<nx; ++i)
            {
                // slightly perturbed ray to avoid hitting edges exactly
                const double xc = ox+(i+0.5)*hx + 1.234567e-7*hx;
                const double yc = oy+(j+0.5)*hy + 2.345678e-7*hy;

                zhit.clear();
                for(const tri& t : T)
                {
                    const double x0=t[0],y0=t[1],x1=t[3],y1=t[4],x2=t[6],y2=t[7];
                    const double den = (y1-y2)*(x0-x2) + (x2-x1)*(y0-y2);
                    if(std::fabs(den)<1.0e-300)
                    continue;
                    const double l0 = ((y1-y2)*(xc-x2) + (x2-x1)*(yc-y2))/den;
                    const double l1 = ((y2-y0)*(xc-x2) + (x0-x2)*(yc-y2))/den;
                    const double l2 = 1.0-l0-l1;
                    if(l0<0.0 || l1<0.0 || l2<0.0)
                    continue;
                    zhit.push_back(l0*t[2] + l1*t[5] + l2*t[8]);
                }

                if(zhit.size()<2)
                continue;

                std::sort(zhit.begin(),zhit.end());

                for(size_t q=0; q+1<zhit.size(); q+=2)
                for(int k=0; k<nz; ++k)
                {
                    const double zc = oz+(k+0.5)*hz;
                    if(zc>zhit[q] && zc<zhit[q+1])
                    vox[voxel(i,j,k)] = mi;
                }
            }
        }
    }
}

void fem_solid::make_nodes_elements()
{
    const int NX=nx+1, NY=ny+1, NZ=nz+1;
    std::vector<int> nid((size_t)NX*NY*NZ,-1);

    auto lattice_node = [&](int i,int j,int k)->int
    {
        const size_t q = (size_t)i + (size_t)NX*((size_t)j + (size_t)NY*(size_t)k);
        if(nid[q]<0)
        {
            nid[q] = (int)X.size();
            X.push_back(Vec3(ox+i*hx,oy+j*hy,oz+k*hz));
        }
        return nid[q];
    };

    X.clear();
    elems.clear();
    vox_elem.assign(vox.size(),-1);

    for(int k=0; k<nz; ++k)
    for(int j=0; j<ny; ++j)
    for(int i=0; i<nx; ++i)
    {
        const int mi = vox[voxel(i,j,k)];
        if(mi<0)
        continue;

        element e;
        e.mat = mi;
        e.ix=i; e.iy=j; e.iz=k;
        e.n[0] = lattice_node(i  ,j  ,k  );
        e.n[1] = lattice_node(i+1,j  ,k  );
        e.n[2] = lattice_node(i+1,j+1,k  );
        e.n[3] = lattice_node(i  ,j+1,k  );
        e.n[4] = lattice_node(i  ,j  ,k+1);
        e.n[5] = lattice_node(i+1,j  ,k+1);
        e.n[6] = lattice_node(i+1,j+1,k+1);
        e.n[7] = lattice_node(i  ,j+1,k+1);

        vox_elem[voxel(i,j,k)] = (int)elems.size();
        elems.push_back(e);
    }

    ngp = full_int ? 8 : 1;
    gps.assign(elems.size()*ngp,gpstate());

    // node -> element connectivity (CSR)
    const int nn = (int)X.size();
    node_elem_start.assign(nn+1,0);
    for(const element& e : elems)
    for(int a=0; a<8; ++a)
    ++node_elem_start[e.n[a]+1];
    for(int i=0; i<nn; ++i)
    node_elem_start[i+1] += node_elem_start[i];
    node_elem.assign(node_elem_start[nn],0);
    std::vector<int> fill(node_elem_start.begin(),node_elem_start.end()-1);
    for(int e=0; e<(int)elems.size(); ++e)
    for(int a=0; a<8; ++a)
    node_elem[fill[elems[e].n[a]]++] = e;
}

void fem_solid::build_surface()
{
    // local faces, ordered so that (x2-x0) x (x3-x1) points outward
    static const int lf[6][4] = {{0,4,7,3},{1,2,6,5},{0,1,5,4},{3,7,6,2},{0,3,2,1},{4,5,6,7}};
    static const int off[6][3] = {{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};

    faces.clear();

    for(int e=0; e<nelem(); ++e)
    {
        const element& el = elems[e];
        if(!el.alive)
        continue;

        for(int f=0; f<6; ++f)
        {
            const int nb = voxel(el.ix+off[f][0],el.iy+off[f][1],el.iz+off[f][2]);
            const int ne = nb<0 ? -1 : vox_elem[nb];
            if(ne>=0 && elems[ne].alive)
            continue;

            face fc;
            for(int q=0; q<4; ++q)
            fc.n[q] = el.n[lf[f][q]];
            fc.elem = e;
            faces.push_back(fc);
        }
    }

    // intact elements per node, debris particles
    nalive.assign(nnode(),0);
    for(const element& el : elems)
    if(el.alive)
    for(int a=0; a<8; ++a)
    ++nalive[el.n[a]];

    orphan.clear();
    for(int i=0; i<nnode(); ++i)
    if(nalive[i]==0 && m[i]>0.0)
    orphan.push_back(i);

    ++surf_version;
    surf_dirty = false;
}

void fem_solid::count_bodies()
{
    // connected components over intact elements (face, edge or corner sharing)
    std::vector<int> comp(nnode(),-1);
    std::vector<int> stack;
    int nc = 0;

    for(int s=0; s<nnode(); ++s)
    {
        if(comp[s]>=0 || nalive[s]==0)
        continue;

        comp[s] = nc;
        stack.push_back(s);

        while(!stack.empty())
        {
            const int i = stack.back();
            stack.pop_back();

            for(int q=node_elem_start[i]; q<node_elem_start[i+1]; ++q)
            {
                const element& el = elems[node_elem[q]];
                if(!el.alive)
                continue;
                for(int a=0; a<8; ++a)
                if(comp[el.n[a]]<0)
                {
                    comp[el.n[a]] = nc;
                    stack.push_back(el.n[a]);
                }
            }
        }
        ++nc;
    }

    nbodies = nc + (int)orphan.size();
}
