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

#include"6DOF_obj.h"
#include"6DOF_obj_remesh.h"
#include"lexer.h"
#include"ghostcell.h"
#include<algorithm>
#include<chrono>

// X 185 4: re-triangulate the 6DOF body surface with an adaptive isotropic remesher.
//
// Target edge length h(x) = X 186 * local cell size of the fluid grid, where the local cell
// size is the smallest spacing of the active directions (x; y in 3D; z unless the grid is a
// single layer). The grid is the global one, assembled from all subdomains, so the result
// is identical for every decomposition. For sigma grids (NHFLOW, FNPF) the vertical node
// positions are mapped to the still water column z = zb + sigma*(wd - zb), zb = global_zmin
// (flat bed assumption under the body).
//
// The surface is remeshed on rank 0 and broadcast: the force integration assigns triangles to
// ranks by their centroid, so all ranks must hold bit-identical triangle lists.
// Sharp edges (X 187 feature angle), corners and open boundaries are preserved; every vertex
// stays on the original STL surface and the triangle orientation is kept.

namespace
{
double grid_cellsize(const std::vector<double> &g, double x)
{
    const int n = int(g.size());

    if(n<2)
    return 1.0e300;

    int i = int(std::upper_bound(g.begin(),g.end(),x) - g.begin()) - 1;
    i = i<0 ? 0 : (i>n-2 ? n-2 : i);

    return g[i+1]-g[i];
}
}

void sixdof_obj::geometry_remesh(lexer *p, ghostcell *pgc)
{
    // ---------------------------------------------------------- global grid node coordinates
    auto gather_nodes = [&](std::vector<double> &g, int gkno, int origin, int kno, const double *XN)
    {
        g.assign(gkno+1,-1.0e300);

        for(int n=0; n<=kno; ++n)
        {
            const int gn = origin + n;

            if(gn>=0 && gn<=gkno)
            g[gn] = XN[n+marge];
        }

        pgc->globalmax(g.data(),gkno+1);
    };

    std::vector<double> gx,gy,gz;

    gather_nodes(gx,p->gknox,p->origin_i,p->knox,p->XN);

    if(p->j_dir==1 && p->gknoy>1)
    gather_nodes(gy,p->gknoy,p->origin_j,p->knoy,p->YN);

    if(p->gknoz>1)
    {
        gather_nodes(gz,p->gknoz,p->origin_k,p->knoz,p->ZN);

        if(p->G2==1)
        {
            const double zb = p->global_zmin;
            const double zt = p->wd>zb ? p->wd : p->global_zmax;

            for(auto &z : gz)
            z = zb + z*(zt-zb);
        }
    }

    const double fac = p->X186;

    auto hfunc = [&](double x, double y, double z)
    {
        double hh = grid_cellsize(gx,x);

        if(!gy.empty())
        hh = std::min(hh,grid_cellsize(gy,y));

        if(!gz.empty())
        hh = std::min(hh,grid_cellsize(gz,z));

        return fac*hh;
    };

    // ---------------------------------------------------------- entity ranges
    int nent = entity_sum;
    bool ranges_ok = (nent>0);
    int pos = 0;

    for(int qn=0; qn<nent && ranges_ok; ++qn)
    {
        if(tstart[qn]!=pos || tend[qn]<tstart[qn])
        ranges_ok = false;

        pos = tend[qn];
    }

    if(pos!=tricount)
    ranges_ok = false;

    if(!ranges_ok)
    {
        // treat all triangles as one surface
        for(int qn=0; qn<nent; ++qn)
        tstart[qn] = tend[qn] = tricount;

        tstart[0] = 0;
        tend[0] = tricount;
        nent = 1;
    }

    // ---------------------------------------------------------- remesh on rank 0
    sixdof_remesh::params prm;
    prm.feature_angle = p->X187;
    prm.iterations = p->X189;

    std::vector<int> ntri(entity_sum,0);
    std::vector<double> buf;

    if(p->mpirank==0)
    {
        auto t0 = std::chrono::steady_clock::now();

        sixdof_remesh R;
        std::vector<sixdof_remesh::vec3> in,out;

        for(int qn=0; qn<nent; ++qn)
        {
            in.clear();

            for(int n=tstart[qn]; n<tend[qn]; ++n)
            for(int q=0; q<3; ++q)
            in.push_back({tri_x[n][q],tri_y[n][q],tri_z[n][q]});

            sixdof_remesh::stats st;
            const bool ok = R.remesh(in,out,hfunc,prm,st);

            cout<<endl<<"6DOF surface remeshing, body "<<n6DOF<<" entity "<<qn<<": "<<(ok?"ok":"FAILED, keeping the input triangles")<<endl;

            if(!ok && st.ntri_estimate>prm.max_tri)
            cout<<"  estimated triangle count "<<st.ntri_estimate<<" exceeds "<<prm.max_tri<<", increase X 186"<<endl;

            sixdof_remesh::print_stats(cout,st);

            if(st.n_boundary_edges>0 || st.n_nonmanifold_edges>0 || st.n_inconsistent_edges>0)
            cout<<"  WARNING: the surface is not closed and consistently oriented, check the STL"<<endl;

            ntri[qn] = int(out.size()/3);

            for(auto &v : out)
            {
                buf.push_back(v[0]);
                buf.push_back(v[1]);
                buf.push_back(v[2]);
            }
        }

        cout<<"  remeshing time: "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count()<<" s"<<endl<<endl;
    }

    // ---------------------------------------------------------- broadcast
    if(entity_sum>0)
    pgc->bcast_int(ntri.data(),entity_sum);

    int total = 0;
    for(int qn=0; qn<nent; ++qn)
    total += ntri[qn];

    buf.resize(9*size_t(total));

    if(total>0)
    pgc->bcast_double(buf.data(),9*total);

    // ---------------------------------------------------------- store
    p->Dresize(tri_x,tricount,total,3,3);
	p->Dresize(tri_y,tricount,total,3,3);
	p->Dresize(tri_z,tricount,total,3,3);
	p->Dresize(tri_x0,tricount,total,3,3);
	p->Dresize(tri_y0,tricount,total,3,3);
	p->Dresize(tri_z0,tricount,total,3,3);

    tricount = total;

    for(int n=0; n<tricount; ++n)
    for(int q=0; q<3; ++q)
    {
        tri_x[n][q] = buf[9*size_t(n)+3*q+0];
        tri_y[n][q] = buf[9*size_t(n)+3*q+1];
        tri_z[n][q] = buf[9*size_t(n)+3*q+2];
    }

    pos = 0;
    for(int qn=0; qn<nent; ++qn)
    {
        tstart[qn] = pos;
        pos += ntri[qn];
        tend[qn] = pos;
    }
    
    for(int qn=nent; qn<entity_sum; ++qn)
    tstart[qn] = tend[qn] = tricount;
}
