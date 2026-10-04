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
Architect: Hans Bihs
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
    {
        stltri[s] = read_stl(stls[s].file);

        // transforms: scale (about the origin), rotation about the vertical axis
        // through the centre of the bounding box, translation
        const shape_stl& st = stls[s];
        if(st.scale!=1.0 || st.rot!=0.0 || st.move.squaredNorm()>0.0)
        {
            double l[3], h[3];
            stl_bounds(stltri[s],l,h);
            const double cx = 0.5*(l[0]+h[0])*st.scale, cy = 0.5*(l[1]+h[1])*st.scale;
            const double ca = std::cos(st.rot*3.14159265358979/180.0), sa = std::sin(st.rot*3.14159265358979/180.0);
            for(tri& t : stltri[s])
            for(int v=0; v<3; ++v)
            {
                double px = t[3*v]*st.scale, py = t[3*v+1]*st.scale, pz = t[3*v+2]*st.scale;
                const double rx = cx + ca*(px-cx) - sa*(py-cy);
                const double ry = cy + sa*(px-cx) + ca*(py-cy);
                t[3*v] = rx + st.move(0);
                t[3*v+1] = ry + st.move(1);
                t[3*v+2] = pz + st.move(2);
            }
        }
    }
    stl_tris.assign(stltri.begin(),stltri.end());

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
    vox_shape.assign((size_t)nx*ny*nz,-1);

    auto matindex = [&](int id)->int
    {
        for(size_t n=0; n<mats.size(); ++n)
        if(mats[n].id==id) return (int)n;
        throw std::runtime_error("FEM: material "+std::to_string(id)+" not defined");
    };

    for(int sc=0; sc<(int)shapes.size(); ++sc)
    {
        const shape_cmd& c = shapes[sc];
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
                {
                    vox[voxel(i,j,k)] = mi;
                    vox_shape[voxel(i,j,k)] = sc;
                }
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
                    {
                        vox[voxel(i,j,k)] = mi;
                        vox_shape[voxel(i,j,k)] = sc;
                    }
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
        e.shape = vox_shape[voxel(i,j,k)];
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

    body = comp;
    body_fixed.assign(nc,0);
    for(int i=0; i<nnode(); ++i)
    if(comp[i]>=0 && fixed[i])
    body_fixed[comp[i]] = 1;
}

// ----------------------------------------------------------------------
// surface snapping: the surface nodes of the voxel mesh are moved onto the
// closest point of the box / STL surface they belong to, which removes the
// stair steps. Elements that would become too distorted keep their nodes.
// ----------------------------------------------------------------------

namespace
{
    // closest point on triangle abc to p (Ericson, Real-Time Collision Detection 5.1.5)
    Eigen::Vector3d closest_on_triangle(const Eigen::Vector3d& p,const Eigen::Vector3d& a,const Eigen::Vector3d& b,const Eigen::Vector3d& c)
    {
        const Eigen::Vector3d ab = b-a, ac = c-a, ap = p-a;
        const double d1 = ab.dot(ap), d2 = ac.dot(ap);
        if(d1<=0.0 && d2<=0.0) return a;
        const Eigen::Vector3d bp = p-b;
        const double d3 = ab.dot(bp), d4 = ac.dot(bp);
        if(d3>=0.0 && d4<=d3) return b;
        const double vc = d1*d4 - d3*d2;
        if(vc<=0.0 && d1>=0.0 && d3<=0.0) return a + (d1/(d1-d3))*ab;
        const Eigen::Vector3d cp = p-c;
        const double d5 = ab.dot(cp), d6 = ac.dot(cp);
        if(d6>=0.0 && d5<=d6) return c;
        const double vb = d5*d2 - d1*d6;
        if(vb<=0.0 && d2>=0.0 && d6<=0.0) return a + (d2/(d2-d6))*ac;
        const double va = d3*d6 - d5*d4;
        if(va<=0.0 && (d4-d3)>=0.0 && (d5-d6)>=0.0) return b + ((d4-d3)/((d4-d3)+(d5-d6)))*(c-b);
        const double denom = 1.0/(va+vb+vc);
        return a + ab*(vb*denom) + ac*(vc*denom);
    }

    Eigen::Vector3d closest_on_box(const Eigen::Vector3d& p,const double* lo,const double* hi)
    {
        Eigen::Vector3d q = p;
        bool inside = true;
        for(int d=0; d<3; ++d)
        {
            if(p(d)<lo[d]) {q(d) = lo[d]; inside = false;}
            if(p(d)>hi[d]) {q(d) = hi[d]; inside = false;}
        }
        if(!inside)
        return q;
        // inside: project onto the nearest face
        int bd = 0; double bdist = 1.0e300, bval = 0.0;
        for(int d=0; d<3; ++d)
        {
            if(p(d)-lo[d]<bdist) {bdist = p(d)-lo[d]; bd = d; bval = lo[d];}
            if(hi[d]-p(d)<bdist) {bdist = hi[d]-p(d); bd = d; bval = hi[d];}
        }
        q(bd) = bval;
        return q;
    }
}

void fem_solid::snap_surface()
{
    static const int lf[6][4] = {{0,4,7,3},{1,2,6,5},{0,1,5,4},{3,7,6,2},{0,3,2,1},{4,5,6,7}};
    static const int off[6][3] = {{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};

    // surface nodes and the shapes whose surface they may belong to
    std::vector<std::vector<int>> cand(nnode());
    for(const element& el : elems)
    for(int f=0; f<6; ++f)
    {
        const int nb = voxel(el.ix+off[f][0],el.iy+off[f][1],el.iz+off[f][2]);
        if(nb>=0 && vox[nb]>=0)
        continue;
        for(int q=0; q<4; ++q)
        {
            std::vector<int>& c = cand[el.n[lf[f][q]]];
            if(el.shape>=0 && std::find(c.begin(),c.end(),el.shape)==c.end()) c.push_back(el.shape);
            // carved by a remove shape
            if(nb>=0 && vox_shape[nb]>=0 && std::find(c.begin(),c.end(),vox_shape[nb])==c.end()) c.push_back(vox_shape[nb]);
        }
    }

    const double hm = std::min(hx,std::min(hy,hz));
    const double dmax = 0.75*hm;

    // triangle bins per STL (cell size 2 h)
    const double cs = 2.0*hm;
    struct bins {double lo[3]; int n[3]; std::vector<std::vector<int>> cell;};
    std::vector<bins> B(stl_tris.size());
    for(size_t s=0; s<stl_tris.size(); ++s)
    {
        const auto& T = stl_tris[s];
        bins& b = B[s];
        double hi[3];
        for(int d=0; d<3; ++d) {b.lo[d] = 1.0e300; hi[d] = -1.0e300;}
        for(const auto& t : T)
        for(int v=0; v<3; ++v)
        for(int d=0; d<3; ++d)
        {
            b.lo[d] = std::min(b.lo[d],t[3*v+d]);
            hi[d] = std::max(hi[d],t[3*v+d]);
        }
        for(int d=0; d<3; ++d)
        {
            b.lo[d] -= dmax;
            b.n[d] = std::max(1,(int)std::ceil((hi[d]+dmax-b.lo[d])/cs));
        }
        b.cell.assign((size_t)b.n[0]*b.n[1]*b.n[2],std::vector<int>());
        for(int k=0; k<(int)T.size(); ++k)
        {
            int i0[3], i1[3];
            for(int d=0; d<3; ++d)
            {
                const double tl = std::min(T[k][d],std::min(T[k][3+d],T[k][6+d])) - dmax;
                const double th = std::max(T[k][d],std::max(T[k][3+d],T[k][6+d])) + dmax;
                i0[d] = std::max(0,(int)std::floor((tl-b.lo[d])/cs));
                i1[d] = std::min(b.n[d]-1,(int)std::floor((th-b.lo[d])/cs));
            }
            for(int a=i0[0]; a<=i1[0]; ++a)
            for(int c=i0[1]; c<=i1[1]; ++c)
            for(int e=i0[2]; e<=i1[2]; ++e)
            b.cell[(size_t)a + (size_t)b.n[0]*((size_t)c + (size_t)b.n[1]*(size_t)e)].push_back(k);
        }
    }

    const std::vector<Vec3> Xlat = X;
    std::vector<char> moved(nnode(),0);

    for(int i=0; i<nnode(); ++i)
    {
        if(cand[i].empty())
        continue;

        const Vec3 p = X[i];
        Vec3 best = p;
        double dbest = 1.0e300;

        for(int sc : cand[i])
        {
            const shape_cmd& c = shapes[sc];
            Vec3 q;
            if(c.type==0)
            {
                const shape_box& bx = boxes[c.idx];
                const double lo[3] = {bx.x0,bx.y0,bx.z0}, hi[3] = {bx.x1,bx.y1,bx.z1};
                q = closest_on_box(p,lo,hi);
            }
            else
            {
                const auto& T = stl_tris[c.idx];
                const bins& b = B[c.idx];
                int ic[3];
                bool out = false;
                for(int d=0; d<3; ++d)
                {
                    ic[d] = (int)std::floor((p(d)-b.lo[d])/cs);
                    if(ic[d]<0 || ic[d]>=b.n[d]) out = true;
                }
                if(out)
                continue;
                double dl = 1.0e300;
                q = p;
                for(int k : b.cell[(size_t)ic[0] + (size_t)b.n[0]*((size_t)ic[1] + (size_t)b.n[1]*(size_t)ic[2])])
                {
                    const Vec3 r = closest_on_triangle(p,Vec3(T[k][0],T[k][1],T[k][2]),Vec3(T[k][3],T[k][4],T[k][5]),Vec3(T[k][6],T[k][7],T[k][8]));
                    const double dd = (r-p).squaredNorm();
                    if(dd<dl) {dl = dd; q = r;}
                }
                if(dl>=1.0e300)
                continue;
            }
            const double dd = (q-p).norm();
            if(dd<dbest) {dbest = dd; best = q;}
        }

        if(dbest<=dmax && dbest>1.0e-12*hm)
        {
            if(plane_strain)
            best(1) = p(1);
            X[i] = best;
            moved[i] = 1;
        }
    }

    // distortion guard: elements must keep a positive Jacobian at all Gauss
    // points and at least 20 % of their volume, otherwise their nodes go back
    const double vvox = hx*hy*hz;
    for(int it=0; it<20; ++it)
    {
        int nbad = 0;
        for(const element& el : elems)
        {
            bool touched = false;
            for(int a=0; a<8; ++a) if(moved[el.n[a]]) touched = true;
            if(!touched)
            continue;
            Vec3 Xa[8];
            for(int a=0; a<8; ++a) Xa[a] = X[el.n[a]];
            egeom G;
            const bool ok = element_geometry(Xa,G);
            double wmin = 1.0e300;
            for(int q=0; q<8; ++q) wmin = std::min(wmin,G.wg[q]);
            if(!ok || G.V<0.2*vvox || wmin<0.05*vvox/8.0)
            {
                for(int a=0; a<8; ++a)
                if(moved[el.n[a]])
                {
                    X[el.n[a]] = Xlat[el.n[a]];
                    moved[el.n[a]] = 0;
                }
                ++nbad;
            }
        }
        if(nbad==0)
        break;
    }
}
