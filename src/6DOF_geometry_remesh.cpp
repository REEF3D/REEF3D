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

#include"6DOF_geometry.h"
#include"6DOF_obj_remesh.h"
#include<cmath>
#include"lexer.h"
#include"ghostcell.h"
#include<algorithm>
#include<chrono>

// X 185 4: re-triangulate the 6DOF body surface with an adaptive isotropic remesher.
//
// Target edge length = X 186 (default 1) * local cell size of the fluid grid, with the cell size
// (X 190) the geometric mean of the active spacings (2, default, regular triangles of one cell
// volume equivalent size), the grid metric (1, an edge spans X 186 cells in every direction,
// triangles stretched like the cells) or the smallest spacing (3). Active directions: x; y in 3D;
// z unless the grid is a single layer. The grid is the global one, assembled from all subdomains, so the result
// is identical for every decomposition. For sigma grids (NHFLOW, FNPF) the vertical node
// positions are mapped to the still water column z = zb + sigma*(wd - zb), zb = global_zmin
// (flat bed assumption under the body).  With mesh refinement (nhflow_amr with G 12, fnpf_amr, the
// CFD mesh refinement with X 10 1) the horizontal spacings are those of the finest level (amr_hfac), with
// the sigma layers doubled on every level (fnpf_amr with G 6 1) also the vertical ones (amr_vfac).
//
// The surface is remeshed on rank 0 and broadcast: the force integration assigns triangles to
// ranks by their centroid, so all ranks must hold bit-identical triangle lists.
// Holes in the STL are closed with X 184 1 (harmonic, flat/membrane patch, default) or X 184 2
// (biharmonic, curvature continuous patch); X 184 0 keeps them as open boundaries.
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

void sixdof_geometry::geometry_remesh(lexer *p, ghostcell *pgc)
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

    // curvilinear grid (G 50): XN/YN are the logical box, the horizontal target length is the mean physical
    // width of the cells (one uniform spacing in x and y)
    if(p->geo_curv==1)
    {
        double ws=0.0;
        int nw=0;
        for(int ii=0; ii<p->knox; ++ii)
        for(int jj=0; jj<p->knoy; ++jj)
        {
            ws += p->geo_width(ii,jj);
            ++nw;
        }
        ws = pgc->globalsum(ws);
        nw = pgc->globalisum(nw);
        // rounded to 6 digits: the remeshing reacts to the last digit of the target length, which the order of
        // the sum over the ranks would change
        double wm = ws/double(std::max(nw,1));
        if(wm>0.0)
        {
            const double sc = std::pow(10.0, 5.0 - std::floor(std::log10(wm)));
            wm = std::round(wm*sc)/sc;
        }
        gx.assign({0.0, wm});
        gy.assign({0.0, wm});
    }

    if(p->gknoz>1)
    {
        gather_nodes(gz,p->gknoz,p->origin_k,p->knoz,p->ZN);

        if(p->G502==1)
        {
            const double zb = p->global_zmin;
            const double zt = p->wd>zb ? p->wd : p->global_zmax;

            for(auto &z : gz)
            z = zb + z*(zt-zb);
        }
    }

    const double fac = p->X186;
    const int mode = p->X190;
    
    // target lengths along x,y,z (one fluid cell times X 186)
    //   X 190 1: grid metric, every edge spans about X 186 cells in each direction (anisotropic)
    //   X 190 2: isotropic, cell size = geometric mean of the active spacings (default)
    //   X 190 3: isotropic, cell size = smallest active spacing
    // inactive directions (2D: y, single layer: z) get the isotropic cell size
    auto Hfunc = [&](double x, double y, double z)
    {
        double d[3] = {amr_hfac*grid_cellsize(gx,x),
                       gy.empty() ? -1.0 : amr_hfac*grid_cellsize(gy,y),
                       gz.empty() ? -1.0 : amr_vfac*grid_cellsize(gz,z)};
        
        double gm=1.0, mn=1.0e300;
        int na=0;
        
        for(int q=0; q<3; ++q)
        if(d[q]>0.0)
        {
            gm *= d[q];
            mn = std::min(mn,d[q]);
            ++na;
        }
        
        gm = std::pow(gm,1.0/double(na));
        
        const double hiso = (mode==3) ? mn : gm;
        
        sixdof_remesh::vec3 H;
        
        for(int q=0; q<3; ++q)
        H[q] = fac*((mode==1 && d[q]>0.0) ? d[q] : hiso);
        
        return H;
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
    prm.hole_fill = p->X184;

    std::vector<int> ntri(entity_sum,0);
    std::vector<double> buf;
    int repaired = 0;

    if(p->mpirank==0)
    {
        auto t0 = std::chrono::steady_clock::now();

        sixdof_remesh R;
        std::vector<sixdof_remesh::vec3> in,out;

        for(int qn=0; qn<nent; ++qn)
        {
            in.clear();

            // curvilinear grid (G 50, map coordinates): remeshed about the centre of the entity's box, so that
            // the surface does not depend on where the body lies (the round-off of coordinates of 1e6 m)
            double ox=0.0, oy=0.0, oz=0.0;
            if(p->geo_curv==1 && tend[qn]>tstart[qn])
            {
                double bx0=1.0e300, bx1=-1.0e300, by0=1.0e300, by1=-1.0e300, bz0=1.0e300, bz1=-1.0e300;
                for(int n=tstart[qn]; n<tend[qn]; ++n)
                for(int q=0; q<3; ++q)
                {
                    bx0 = std::min(bx0,tri_x[n][q]); bx1 = std::max(bx1,tri_x[n][q]);
                    by0 = std::min(by0,tri_y[n][q]); by1 = std::max(by1,tri_y[n][q]);
                    bz0 = std::min(bz0,tri_z[n][q]); bz1 = std::max(bz1,tri_z[n][q]);
                }
                ox = 0.5*(bx0+bx1);
                oy = 0.5*(by0+by1);
                oz = 0.5*(bz0+bz1);
            }

            for(int n=tstart[qn]; n<tend[qn]; ++n)
            for(int q=0; q<3; ++q)
            in.push_back({tri_x[n][q]-ox,tri_y[n][q]-oy,tri_z[n][q]-oz});

            auto Hloc = [&](double x, double y, double z) { return Hfunc(x+ox, y+oy, z+oz); };

            sixdof_remesh::stats st;
            const bool ok = R.remesh(in,out,sixdof_remesh::metric_func(Hloc),prm,st);

            for(auto &v : out)
            {
                v[0] += ox;
                v[1] += oy;
                v[2] += oz;
            }

            cout<<endl<<"6DOF surface remeshing, body "<<id<<" entity "<<qn<<": "<<(ok?"ok":"FAILED, keeping the input triangles")<<endl;

            if(!ok && st.ntri_estimate>prm.max_tri)
            cout<<"  estimated triangle count "<<st.ntri_estimate<<" exceeds "<<prm.max_tri<<", increase X 186"<<endl;

            sixdof_remesh::print_stats(cout,st);

            if(st.n_holes_filled>0 || st.n_gap_merges>0)
            {
                cout<<"  NOTE: the STL surface had holes or gaps, they were closed (X 184 "<<p->X184<<")"<<endl;
                repaired = 1;
            }
            
            if(st.n_boundary_edges>0 || st.n_nonmanifold_edges>0 || st.n_inconsistent_edges>0)
            cout<<"  WARNING: the surface is not closed and consistently oriented, check the STL"
                <<(p->X184==0 ? " or activate hole filling with X 184 1" : "")<<endl;

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
    
    pgc->bcast_int(&repaired,1);

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
    
    // the volume of an STL with holes is meaningless: the owner recomputes volume and mass (X 21)
    // or density (X 22) from the closed surface
    holes_closed = (repaired==1);
}
