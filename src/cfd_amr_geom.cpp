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

#include"cfd_amr.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"geo_mesh.h"
#include"geo_raycast.h"
#include"geo_grid_solids.h"
#include<cmath>
#include<iostream>
#include<limits>

//  The solids and the topography of the grid file (S, T entities of DIVEMesh, the geometry core
//  geo_mesh) on the patches, at the resolution of the patch.
//
//  Level 0 casts the entities on its grid (lexer::grid_solids: inside/outside by ray parity, the
//  distance to the crossings along the grid lines) and reinitialises the field (reinitopo).  A
//  new patch first takes the field of its parent (prolongation, the registry: solid, topo), whose
//  zero level lies where level 0 has put it.  The entities are then cast on the patch: the sign of
//  every cell from the ray parity on the patch grid (lexer::grid_solids_io), the distance within a
//  band of geo_band cells of the patch from the triangles themselves (the exact point-triangle
//  distance; in 2D the distance to the cut of the triangles with the plane of the grid), outside
//  the band the prolonged distance (at least the band width).  The direct forcing of the patch
//  (solid_forcing_ini: Hsolidface, fbh5 within 1.1 dx) and the flags G 25 read the field near the
//  surface, so the patch sees the geometry at its own resolution: a cylinder of 4 level-0 cells
//  is a cylinder on level 2, not a polygon of level 0.
//
//  As on level 0 a face of a solid on the domain boundary is no fluid-solid interface (the
//  triangles in a plane of the boundary are left out of the distance of the solids).
//
//  The geodat bed (G 9, G 10 of DIVEMesh: the bed level of every level-0 column, GEOB) is a surface
//  through the centres of the level-0 columns, linear on two triangles per four neighbouring
//  centres (in 2D a polyline), held level from the outermost centres to the domain boundary and
//  beyond (as level 0 takes the edge column there).  A cell is in the bed where it lies below that
//  surface; its distance within the band is the exact distance to the triangles.  Level 0 has the
//  vertical distance to the bed level of its column, reinitialised, so on a slope or at a step of
//  the bed the patch is no longer a copy of that.  The bed has no finer data than the level-0
//  columns: the surface is the one the level-0 columns define, not a finer bathymetry.  The levels of
//  all columns are gathered once (gbed0), so a patch reads the columns around it on any rank (G 40).
//  With the entities of the same role the field is the union, the minimum of both, as on level 0.
//  S 57 (a flat initial bed of the topography) keeps the prolonged field.
//
//  The cells around the patch keep what the fills give them (a patch of the same level, or the
//  coarser level), the cells outside the domain the boundary values (bc_state).

namespace
{
const int geo_band = 5;         // cells of the patch within which the exact distance is taken

// [s,e): the cells of [lo,hi) whose centre lies in [x0,x1]
void geo_range(const double *XP, int lo, int hi, double x0, double x1, int &s, int &e)
{
    const int m = increment::marge;
    int a = lo, b = hi;
    while(a<b)
    {
        const int c = (a+b)/2;
        if(XP[c+m]<x0)
        a = c+1;
        else
        b = c;
    }
    s = a;
    a = s; b = hi;
    while(a<b)
    {
        const int c = (a+b)/2;
        if(XP[c+m]<=x1)
        a = c+1;
        else
        b = c;
    }
    e = a;
}

// squared distance of (px,pz) to the segment (ax,az)-(bx,bz)
double dist2_seg(double px, double pz, double ax, double az, double bx, double bz)
{
    const double ux = bx-ax, uz = bz-az;
    const double L2 = ux*ux + uz*uz;
    double t = (L2>0.0) ? ((px-ax)*ux + (pz-az)*uz)/L2 : 0.0;
    t = (t<0.0) ? 0.0 : (t>1.0 ? 1.0 : t);
    const double dx = px-(ax+t*ux), dz = pz-(az+t*uz);
    return dx*dx + dz*dz;
}
}

// the roles cast on the patches (the constructor): the entities of a role read by the run, its geodat bed
void cfd_amr::geom_setup(lexer *p)
{
    geo_role[0] = geo_role[1] = 0;
    geo_bed[0] = geo_bed[1] = 0;
    if(!solids || p->gridgeo==nullptr)
    return;

    const geo_mesh &G = *p->gridgeo;
    const int roles[2] = {geo_mesh::role_solid, geo_mesh::role_topo};
    const int read[2] = {p->solidread, p->toporead};
    for(int r=0; r<2; ++r)
    {
        bool objects = false;
        for(const geo_object &ob : G.obj)
        if(ob.role==roles[r] && ob.te>ob.ts)
        objects = true;

        const bool bed = (G.geodat==roles[r]);

        if(read[r]==0 || (!objects && !bed))
        continue;

        // S 57: a flat initial bed replaces the topography of the grid file on level 0 (geotopo::dat)
        if(r==1 && p->S57>-1.0e20)
        {
            if(p->mpirank==0)
            cout<<"CFD AMR: the topography with S 57: the patches keep the prolonged distance"<<endl;
            continue;
        }

        geo_role[r] = 1;
        geo_bed[r] = bed ? 1 : 0;
    }

    // the geodat bed of all level-0 columns (every column on one rank: no partition in z)
    if(geo_bed[0]==1 || geo_bed[1]==1)
    {
        const size_t n = (size_t)p->gknox*p->gknoy;
        vector<double> v(n,0.0);
        for(int i=0; i<p->knox; ++i)
        for(int j=0; j<p->knoy; ++j)
        v[(size_t)(i+p->origin_i)*p->gknoy + (j+p->origin_j)] = p->geobed[(i-p->imin)*p->jmax + (j-p->jmin)];
        gbed0.assign(n,0.0);
        MPI_Allreduce(v.data(),gbed0.data(),(int)n,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    }
}

// the entities and the geodat bed on patch id at its resolution (after its state is prolonged or
// copied, before the fill of its level)
void cfd_amr::patch_geometry(int id)
{
    if((geo_role[0]==0 && geo_role[1]==0) || is_shadow(id))
    return;

    cfd_amr_patch *c = CP(id);
    lexer *pp = c->pp;
    fdm *fa = c->a;
    const geo_mesh &G = *pp->gridgeo;
    const int m = marge;
    const int nf = pp->imax*pp->jmax*pp->kmax;
    const bool d2 = (pp->j_dir==0);

    // the interior of the patch
    geo_cart g = geo_raycast::extended(pp,(G.dxm>0.0) ? G.dxm : p0->DXM);
    g.is = 0; g.ie = pp->knox;
    g.js = 0; g.je = pp->knoy;
    g.ks = 0; g.ke = pp->knoz;

    // the band: geo_band cells of the patch
    double h = 0.0;
    for(int i=0; i<pp->knox; ++i)
    h = MAX(h,pp->DXN[i+m]);
    for(int k=0; k<pp->knoz; ++k)
    h = MAX(h,pp->DZN[k+m]);
    if(!d2)
    for(int j=0; j<pp->knoy; ++j)
    h = MAX(h,pp->DYN[j+m]);
    const double band = double(geo_band)*h;
    const double band2 = band*band;

    auto cq = [&](int i, int j, int k) { return (i-pp->imin)*pp->jmax*pp->kmax + (j-pp->jmin)*pp->kmax + k-pp->kmin; };

    // a triangle in a plane of the domain boundary
    const double eps = 1.0e-6*h;
    auto on_boundary = [&](const double *tx, const double *ty, const double *tz)
    {
        auto plane = [&](const double *t, double b) { return fabs(t[0]-b)<eps && fabs(t[1]-b)<eps && fabs(t[2]-b)<eps; };
        return plane(tx,p0->global_xmin) || plane(tx,p0->global_xmax) || plane(ty,p0->global_ymin) || plane(ty,p0->global_ymax)
            || plane(tz,p0->global_zmin) || plane(tz,p0->global_zmax);
    };

    const int roles[2] = {geo_mesh::role_solid, geo_mesh::role_topo};
    vector<int> IO;
    vector<double> dmin, phi;
    const double big = 1.0e30;     // no surface within the band

    for(int r=0; r<2; ++r)
    {
        if(geo_role[r]==0)
        continue;

        const bool obj = geo_grid_solids_io(pp,roles[r],g,IO);
        if(!obj && geo_bed[r]==0)
        continue;

        // the signed distance within the band of the patch cells (+-big outside it): the entities,
        // united with the bed
        phi.assign(nf,big);

        if(obj)
        {
        dmin.assign(nf,std::numeric_limits<double>::max());

        // the triangles of the role: the cells within the band of each
        for(const geo_object &ob : G.obj)
        if(ob.role==roles[r] && ob.te>ob.ts)
        for(int n=ob.ts; n<ob.te; ++n)
        {
            const double *tx = G.tri_x[n], *ty = G.tri_y[n], *tz = G.tri_z[n];

            // solids: a face on the domain boundary is no fluid-solid interface (as grid_solids)
            if(r==0 && on_boundary(tx,ty,tz))
            continue;
            const double x0 = MIN3(tx[0],tx[1],tx[2])-band, x1 = MAX3(tx[0],tx[1],tx[2])+band;
            const double z0 = MIN3(tz[0],tz[1],tz[2])-band, z1 = MAX3(tz[0],tz[1],tz[2])+band;
            int is,ie,js,je,ks,ke;
            geo_range(pp->XP,0,pp->knox,x0,x1,is,ie);
            geo_range(pp->ZP,0,pp->knoz,z0,z1,ks,ke);
            if(ie<=is || ke<=ks)
            continue;

            if(d2)
            {
                // the cut of the triangle with the plane y = YP of the grid: a segment in x-z
                const double y = pp->YP[0+m];
                double sx[3], sz[3];
                int ns = 0;
                for(int e=0; e<3 && ns<3; ++e)
                {
                    const int a = e, b = (e+1)%3;
                    const double fa_ = ty[a]-y, fb_ = ty[b]-y;
                    if(fa_==0.0)
                    {
                        sx[ns] = tx[a]; sz[ns] = tz[a]; ++ns;
                    }
                    else if(fa_*fb_<0.0)
                    {
                        const double t = fa_/(fa_-fb_);
                        sx[ns] = tx[a]+t*(tx[b]-tx[a]); sz[ns] = tz[a]+t*(tz[b]-tz[a]); ++ns;
                    }
                }
                if(ns<2)
                continue;
                for(int i=is; i<ie; ++i)
                for(int k=ks; k<ke; ++k)
                {
                    const double px = pp->XP[i+m], pz = pp->ZP[k+m];
                    double dd = dist2_seg(px,pz,sx[0],sz[0],sx[1],sz[1]);
                    if(ns==3)
                    {
                        dd = MIN(dd,dist2_seg(px,pz,sx[1],sz[1],sx[2],sz[2]));
                        dd = MIN(dd,dist2_seg(px,pz,sx[0],sz[0],sx[2],sz[2]));
                    }
                    for(int j=0; j<pp->knoy; ++j)
                    {
                        double &v = dmin[cq(i,j,k)];
                        v = MIN(v,dd);
                    }
                }
                continue;
            }

            const double y0 = MIN3(ty[0],ty[1],ty[2])-band, y1 = MAX3(ty[0],ty[1],ty[2])+band;
            geo_range(pp->YP,0,pp->knoy,y0,y1,js,je);
            if(je<=js)
            continue;
            for(int i=is; i<ie; ++i)
            for(int j=js; j<je; ++j)
            for(int k=ks; k<ke; ++k)
            {
                const double dd = geo_raycast::dist2_tri(pp->XP[i+m],pp->YP[j+m],pp->ZP[k+m],
                                                         tx[0],ty[0],tz[0],tx[1],ty[1],tz[1],tx[2],ty[2],tz[2]);
                double &v = dmin[cq(i,j,k)];
                v = MIN(v,dd);
            }
        }

        // the sign of the ray parity, the exact distance in the band
        for(int i=0; i<pp->knox; ++i)
        for(int j=0; j<pp->knoy; ++j)
        for(int k=0; k<pp->knoz; ++k)
        {
            const int q = cq(i,j,k);
            phi[q] = ((IO[q]==-1) ? -1.0 : 1.0)*((dmin[q]<band2) ? sqrt(dmin[q]) : big);
        }
        }

        // the geodat bed
        if(geo_bed[r]==1)
        patch_bed(c,band,phi);

        // the field: the distance in the band, the prolonged one (at least the band width) outside it
        field &f = (r==0) ? fa->solid : fa->topo;
        for(int i=0; i<pp->knox; ++i)
        for(int j=0; j<pp->knoy; ++j)
        for(int k=0; k<pp->knoz; ++k)
        {
            const int q = cq(i,j,k);
            if(fabs(phi[q])<0.5*big)
            f.V[q] = phi[q];
            else
            f.V[q] = ((phi[q]<0.0) ? -1.0 : 1.0)*MAX(fabs(f.V[q]),band);
        }
    }
}

// the geodat bed on patch c: the surface through the centres of the level-0 columns (two triangles per
// four centres; 2D: a polyline), held level beyond the outermost centres.  Every interior cell: below
// the surface inside; within the band the exact distance to it; phi = min(phi, that) (the union)
void cfd_amr::patch_bed(cfd_amr_patch *c, double band, vector<double> &phi)
{
    lexer *pp = c->pp;
    const int m = marge;
    const bool d2 = (pp->j_dir==0);
    const double band2 = band*band;
    const double big = 1.0e30;
    const int gny = p0->gknoy;

    auto cq = [&](int i, int j, int k) { return (i-pp->imin)*pp->jmax*pp->kmax + (j-pp->jmin)*pp->kmax + k-pp->kmin; };

    // the vertices along direction d covering [lo,hi]: the level-0 column centres there, one more on
    // each side, and a vertex far outside where the domain ends (the edge column held)
    auto axis = [&](int d, double lo, double hi, vector<double> &vx, vector<int> &vi)
    {
        const int n = gn0[d];
        auto ctr = [&](int I) { return 0.5*(gx[d][I+marge]+gx[d][I+1+marge]); };
        const double far = (hi-lo) + (gx[d][n+marge]-gx[d][marge]) + 1.0;
        vx.clear();
        vi.clear();
        int I0 = 0;
        while(I0<n-1 && ctr(I0+1)<=lo)
        ++I0;
        if(ctr(I0)>lo)
        {
            vx.push_back(MIN(lo,ctr(0))-far);
            vi.push_back(0);
        }
        int I = I0;
        for(; I<n; ++I)
        {
            vx.push_back(ctr(I));
            vi.push_back(I);
            if(ctr(I)>=hi)
            break;
        }
        if(I==n)
        {
            vx.push_back(MAX(hi,ctr(n-1))+far);
            vi.push_back(n-1);
        }
    };

    vector<double> vx, vy;
    vector<int> vi, vj;
    axis(0,pp->XN[m]-band,pp->XN[pp->knox+m]+band,vx,vi);
    if(d2)
    {
        vy.assign(1,pp->YP[0+m]);
        vj.assign(1,0);
    }
    else
    axis(1,pp->YN[m]-band,pp->YN[pp->knoy+m]+band,vy,vj);

    const int nx = (int)vx.size(), ny = (int)vy.size();
    auto zb = [&](int a, int b) { return gbed0[(size_t)vi[a]*gny + vj[b]]; };

    // the interval of v holding x (the vertices cover the patch)
    auto locate = [](const vector<double> &v, double x)
    {
        int a = 0, b = (int)v.size()-1;
        while(b-a>1)
        {
            const int c = (a+b)/2;
            if(v[c]<=x)
            a = c;
            else
            b = c;
        }
        return a;
    };

    // the exact distance of the cells within the band of the surface
    vector<double> d2min((size_t)pp->imax*pp->jmax*pp->kmax,std::numeric_limits<double>::max());

    if(d2)
    {
        for(int a=0; a<nx-1; ++a)
        {
            const double ax = vx[a], az = zb(a,0), bx = vx[a+1], bz = zb(a+1,0);
            int is,ie,ks,ke;
            geo_range(pp->XP,0,pp->knox,MIN(ax,bx)-band,MAX(ax,bx)+band,is,ie);
            geo_range(pp->ZP,0,pp->knoz,MIN(az,bz)-band,MAX(az,bz)+band,ks,ke);
            for(int i=is; i<ie; ++i)
            for(int k=ks; k<ke; ++k)
            {
                const double dd = dist2_seg(pp->XP[i+m],pp->ZP[k+m],ax,az,bx,bz);
                for(int j=0; j<pp->knoy; ++j)
                {
                    double &v = d2min[cq(i,j,k)];
                    v = MIN(v,dd);
                }
            }
        }
    }
    else
    {
        for(int a=0; a<nx-1; ++a)
        for(int b=0; b<ny-1; ++b)
        {
            // (a,b), (a+1,b), (a+1,b+1) and (a,b), (a+1,b+1), (a,b+1)
            const double X[4] = {vx[a],vx[a+1],vx[a+1],vx[a]};
            const double Y[4] = {vy[b],vy[b],vy[b+1],vy[b+1]};
            const double Z[4] = {zb(a,b),zb(a+1,b),zb(a+1,b+1),zb(a,b+1)};
            const int T[2][3] = {{0,1,2},{0,2,3}};
            for(int t=0; t<2; ++t)
            {
                const int *v = T[t];
                int is,ie,js,je,ks,ke;
                geo_range(pp->XP,0,pp->knox,MIN3(X[v[0]],X[v[1]],X[v[2]])-band,MAX3(X[v[0]],X[v[1]],X[v[2]])+band,is,ie);
                geo_range(pp->YP,0,pp->knoy,MIN3(Y[v[0]],Y[v[1]],Y[v[2]])-band,MAX3(Y[v[0]],Y[v[1]],Y[v[2]])+band,js,je);
                geo_range(pp->ZP,0,pp->knoz,MIN3(Z[v[0]],Z[v[1]],Z[v[2]])-band,MAX3(Z[v[0]],Z[v[1]],Z[v[2]])+band,ks,ke);
                for(int i=is; i<ie; ++i)
                for(int j=js; j<je; ++j)
                for(int k=ks; k<ke; ++k)
                {
                    const double dd = geo_raycast::dist2_tri(pp->XP[i+m],pp->YP[j+m],pp->ZP[k+m],
                                                             X[v[0]],Y[v[0]],Z[v[0]],X[v[1]],Y[v[1]],Z[v[1]],X[v[2]],Y[v[2]],Z[v[2]]);
                    double &w = d2min[cq(i,j,k)];
                    w = MIN(w,dd);
                }
            }
        }
    }

    // the sign: below the surface (on the triangle of the cell, the polyline in 2D)
    for(int i=0; i<pp->knox; ++i)
    for(int j=0; j<pp->knoy; ++j)
    {
        const double x = pp->XP[i+m], y = pp->YP[j+m];
        const int a = locate(vx,x);
        double z;
        if(d2)
        {
            const double u = (x-vx[a])/(vx[a+1]-vx[a]);
            z = zb(a,0) + u*(zb(a+1,0)-zb(a,0));
        }
        else
        {
            const int b = locate(vy,y);
            const double u = (x-vx[a])/(vx[a+1]-vx[a]);
            const double w = (y-vy[b])/(vy[b+1]-vy[b]);
            if(u>=w)
            z = zb(a,b) + u*(zb(a+1,b)-zb(a,b)) + w*(zb(a+1,b+1)-zb(a+1,b));
            else
            z = zb(a,b) + w*(zb(a,b+1)-zb(a,b)) + u*(zb(a+1,b+1)-zb(a,b+1));
        }

        for(int k=0; k<pp->knoz; ++k)
        {
            const int q = cq(i,j,k);
            const double s = (pp->ZP[k+m]<z) ? -1.0 : 1.0;
            const double pb = s*((d2min[q]<band2) ? sqrt(d2min[q]) : big);
            phi[q] = MIN(phi[q],pb);
        }
    }
}
