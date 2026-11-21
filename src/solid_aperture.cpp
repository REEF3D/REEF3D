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

#include"solid_aperture.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"increment.h"
#include"solid_closed.h"
#include"field.h"
#include"geo_mesh.h"
#include"geo_raycast.h"
#include"geo_grid_solids.h"
#include<vector>
#include<limits>
#include<algorithm>
#include<cmath>
#include<cstring>

namespace
{
const int solid_band = 3;       // cells within which the exact distance to the solids is taken
const int press_layers = 4;     // layers of closed cells that take the pressure of the fluid

inline double clip(double A)
{
    return (A<solid_aperture_min) ? 0.0 : A;
}

// [s,e): the cells of [lo,hi) whose centre lies in [x0,x1]
void cell_range(const double *XP, int lo, int hi, double x0, double x1, int &s, int &e)
{
    const int m = increment::marge;
    int a = lo, b = hi;
    while(a<b)
    {
        const int c = a + (b-a)/2;
        if(XP[c+m]<x0)
        a = c+1;
        else
        b = c;
    }
    s = a;
    b = hi;
    while(a<b)
    {
        const int c = a + (b-a)/2;
        if(XP[c+m]<=x1)
        a = c+1;
        else
        b = c;
    }
    e = a;
}

// squared distance of (px,pz) to the segment (ax,az)-(bx,bz)
double seg_dist2(double px, double pz, double ax, double az, double bx, double bz)
{
    const double ex = bx-ax, ez = bz-az;
    const double l2 = ex*ex + ez*ez;
    double t = (l2>0.0) ? ((px-ax)*ex + (pz-az)*ez)/l2 : 0.0;
    t = std::max(0.0,std::min(1.0,t));
    const double dx = px-(ax+t*ex), dz = pz-(az+t*ez);
    return dx*dx + dz*dz;
}

// The signed distance of the cell centres (ghost cells as well) to the solids of the grid file
// (S, the triangles of the geometry core) within solid_band cells of their surface, +-the band
// beyond it; false: no solids there.
//
// The solid level set of the grid (a->solid, lexer::grid_solids) is the distance along the grid
// lines to the crossings of the surface, inside a solid not below half a cell: next to an edge of
// a solid and inside it the level set is no distance, and the corners of the faces interpolated
// from it are off (a wall at 30 % of a cell: 46 % open).  The sign from the ray parity of the cell
// centre, the distance to the triangles (2D: to their cut with the plane of the grid), as the
// patches of the mesh refinement take it (cfd_amr::patch_geometry)
bool solid_distance(lexer *p, std::vector<double> &sd)
{
    if(p->gridgeo==nullptr)
    return false;

    const geo_mesh &G = *p->gridgeo;
    const int m = increment::marge;
    const int nf = p->imax*p->jmax*p->kmax;
    const bool d2 = (p->j_dir==0);

    geo_cart g = geo_raycast::extended(p,(G.dxm>0.0) ? G.dxm : p->DXM);

    std::vector<int> IO;
    if(!geo_grid_solids_io(p,geo_mesh::role_solid,g,IO))
    return false;

    const int i0 = p->imin, i1 = p->imin+p->imax;
    const int j0 = d2 ? 0 : p->jmin, j1 = d2 ? 1 : p->jmin+p->jmax;
    const int k0 = p->kmin, k1 = p->kmin+p->kmax;

    double h = 0.0;
    for(int i=0; i<p->knox; ++i)
    h = std::max(h,p->DXN[i+m]);
    for(int k=0; k<p->knoz; ++k)
    h = std::max(h,p->DZN[k+m]);
    if(!d2)
    for(int j=0; j<p->knoy; ++j)
    h = std::max(h,p->DYN[j+m]);

    const double band = double(solid_band)*h;
    const double band2 = band*band;

    auto cq = [&](int i, int j, int k) { return (i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin; };

    // a triangle in a plane of the domain boundary is no surface between fluid and solid
    const double eps = 1.0e-6*h;
    auto on_boundary = [&](const double *tx, const double *ty, const double *tz)
    {
        auto plane = [&](const double *t, double b) { return fabs(t[0]-b)<eps && fabs(t[1]-b)<eps && fabs(t[2]-b)<eps; };
        return plane(tx,p->global_xmin) || plane(tx,p->global_xmax) || (!d2 && (plane(ty,p->global_ymin) || plane(ty,p->global_ymax)))
            || plane(tz,p->global_zmin) || plane(tz,p->global_zmax);
    };

    std::vector<double> dmin(nf,std::numeric_limits<double>::max());

    for(const geo_object &ob : G.obj)
    if(ob.role==geo_mesh::role_solid && ob.te>ob.ts)
    for(int n=ob.ts; n<ob.te; ++n)
    {
        const double *tx = G.tri_x[n], *ty = G.tri_y[n], *tz = G.tri_z[n];

        if(on_boundary(tx,ty,tz))
        continue;

        int is,ie,js,je,ks,ke;
        cell_range(p->XP,i0,i1,MIN3(tx[0],tx[1],tx[2])-band,MAX3(tx[0],tx[1],tx[2])+band,is,ie);
        cell_range(p->ZP,k0,k1,MIN3(tz[0],tz[1],tz[2])-band,MAX3(tz[0],tz[1],tz[2])+band,ks,ke);
        if(ie<=is || ke<=ks)
        continue;

        if(d2)
        {
            // the cut of the triangle with the plane y = YP of the grid: a segment in x-z
            const double y = p->YP[0+m];
            double sx[3], sz[3];
            int ns = 0;
            for(int e=0; e<3 && ns<3; ++e)
            {
                const int a = e, b = (e+1)%3;
                const double fa = ty[a]-y, fb = ty[b]-y;
                if(fa==0.0)
                {
                    sx[ns] = tx[a]; sz[ns] = tz[a]; ++ns;
                }
                else if(fa*fb<0.0)
                {
                    const double t = fa/(fa-fb);
                    sx[ns] = tx[a]+t*(tx[b]-tx[a]); sz[ns] = tz[a]+t*(tz[b]-tz[a]); ++ns;
                }
            }
            if(ns<2)
            continue;

            for(int i=is; i<ie; ++i)
            for(int k=ks; k<ke; ++k)
            {
                const double px = p->XP[i+m], pz = p->ZP[k+m];
                double dd = seg_dist2(px,pz,sx[0],sz[0],sx[1],sz[1]);
                if(ns==3)
                {
                    dd = std::min(dd,seg_dist2(px,pz,sx[1],sz[1],sx[2],sz[2]));
                    dd = std::min(dd,seg_dist2(px,pz,sx[0],sz[0],sx[2],sz[2]));
                }
                double &v = dmin[cq(i,0,k)];
                v = std::min(v,dd);
            }
            continue;
        }

        cell_range(p->YP,j0,j1,MIN3(ty[0],ty[1],ty[2])-band,MAX3(ty[0],ty[1],ty[2])+band,js,je);
        if(je<=js)
        continue;

        for(int i=is; i<ie; ++i)
        for(int j=js; j<je; ++j)
        for(int k=ks; k<ke; ++k)
        {
            const double dd = geo_raycast::dist2_tri(p->XP[i+m],p->YP[j+m],p->ZP[k+m],
                                                     tx[0],ty[0],tz[0],tx[1],ty[1],tz[1],tx[2],ty[2],tz[2]);
            double &v = dmin[cq(i,j,k)];
            v = std::min(v,dd);
        }
    }

    sd.assign(nf,band);
    for(int i=i0; i<i1; ++i)
    for(int j=j0; j<j1; ++j)
    for(int k=k0; k<k1; ++k)
    {
        const int q = cq(i,j,k);
        sd[q] = ((IO[q]==-1) ? -1.0 : 1.0)*((dmin[q]<band2) ? sqrt(dmin[q]) : band);
    }

    return true;
}
}

void solid_aperture_update(lexer *p, fdm *a, ghostcell *pgc)
{
    if(p->D1!=2 || a->apx.V==nullptr)
    return;

    const bool solids = p->solidread>0;
    const bool bed = p->topoforcing>0;

    int i,j,k;
    const int marge = increment::marge;
    const int ncell = p->imax*p->jmax*p->kmax;

    // changed since the last call?  (bitwise: a moving bed changes the sum)
    unsigned long long sum = 0;
    {
        unsigned long long w;
        for(int q=0; q<ncell; ++q)
        {
            if(solids)
            {
            std::memcpy(&w,&a->solid.V[q],sizeof(double));
            sum = sum*1099511628211ULL + w;
            }
            if(bed)
            {
            std::memcpy(&w,&a->topo.V[q],sizeof(double));
            sum = sum*1099511628211ULL + w;
            }
        }
    }

    // (all ranks together: the halo update below is collective.  A patch of the mesh refinement is
    // rank local)
    int changed = (a->ap_done && a->ap_sum==sum) ? 0 : 1;
    if(p->amr_patch==0)
    changed = pgc->globalimax(changed);
    if(changed==0)
    return;
    a->ap_done = true;
    a->ap_sum = sum;

    std::fill(a->apx.V, a->apx.V + ncell, 1.0);
    std::fill(a->apy.V, a->apy.V + ncell, 1.0);
    std::fill(a->apz.V, a->apz.V + ncell, 1.0);

    if(!solids && !bed)
    return;

    const bool jd = p->j_dir==1;

    // the solids: the exact distance to their surface (the level set of the grid without the
    // geometry of the grid file)
    std::vector<double> sd;
    const bool exact = solids && solid_distance(p,sd);

    // the cells of the domain (the ghost cells beyond its edges hold the boundary condition of the
    // level set: the cell next to the edge instead), within the arrays
    const int i0 = p->imin, i1 = p->imin+p->imax-1;
    const int j0 = p->jmin, j1 = p->jmin+p->jmax-1;
    const int k0 = p->kmin, k1 = p->kmin+p->kmax-1;
    const int ilo = (p->periodic1==1) ? i0 : std::max(i0,-p->origin_i), ihi = (p->periodic1==1) ? i1 : std::min(i1,p->gknox-1-p->origin_i);
    const int jlo = (p->periodic2==1) ? j0 : std::max(j0,-p->origin_j), jhi = (p->periodic2==1) ? j1 : std::min(j1,p->gknoy-1-p->origin_j);
    const int klo = (p->periodic3==1) ? k0 : std::max(k0,-p->origin_k), khi = (p->periodic3==1) ? k1 : std::min(k1,p->gknoz-1-p->origin_k);

    // the level set of the solid and the bed, slightly lowered: a surface on a cell face is closed
    const double eps = 1.0e-8*p->DXM;

    auto phi = [&](int ii, int jj, int kk)
    {
        ii = std::max(ilo,std::min(ihi,ii));
        kk = std::max(klo,std::min(khi,kk));
        if(jd)
        jj = std::max(jlo,std::min(jhi,jj));
        else
        jj = 0;

        const double fs = !solids ? 0.0 : (exact ? sd[(ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax + kk-p->kmin] : a->solid(ii,jj,kk));

        double f;
        if(solids && bed)
        f = std::min(fs,a->topo(ii,jj,kk));
        else if(solids)
        f = fs;
        else
        f = a->topo(ii,jj,kk);

        return f - eps;
    };

    // the cell widths, within the arrays
    auto dxn = [&](int ii) { return p->DXN[std::max(i0,std::min(i1,ii))+marge]; };
    auto dyn = [&](int jj) { return p->DYN[std::max(j0,std::min(j1,jj))+marge]; };
    auto dzn = [&](int kk) { return p->DZN[std::max(k0,std::min(k1,kk))+marge]; };

    // interpolation weights of the cell and its upper neighbour to the face between them
    auto wx = [&](int ii) { return dxn(ii+1)/(dxn(ii)+dxn(ii+1)); };
    auto wy = [&](int jj) { return jd ? dyn(jj+1)/(dyn(jj)+dyn(jj+1)) : 1.0; };
    auto wz = [&](int kk) { return dzn(kk+1)/(dzn(kk)+dzn(kk+1)); };

    // the level set at the corner (ii+1/2, jj+1/2, kk+1/2); 2D: at (ii+1/2, kk+1/2), interpolated
    // from the cells around it.  At an edge of the solid the level set has a kink and the
    // interpolation is too high (the corner of a box: the cell diagonally outside adds its distance,
    // the corner on the surface is found in the fluid, and the two faces of the box cell at it open
    // completely: a flow around the edge through the solid).  The means of the cells opposite each
    // other across the corner are equal where the level set is smooth (a curved surface of radius R:
    // they differ by h^2/(2R) at most).  Where they differ by more than solid_kink of the cell width,
    // the corner takes the lowest of them, lowered by solid_kink_low times the difference: towards
    // the solid (the corner of a box from DIVEMesh, S 10: lowest mean 0.26 h, difference 0.24 h)
    const double solid_kink = 0.15;
    const double solid_kink_low = 1.5;

    auto corner = [&](int ii, int jj, int kk)
    {
        const double ax = wx(ii), az = wz(kk);

        if(!jd)
        {
        const double f00 = phi(ii,jj,kk), f10 = phi(ii+1,jj,kk), f01 = phi(ii,jj,kk+1), f11 = phi(ii+1,jj,kk+1);
        const double m1 = 0.5*(f00+f11), m2 = 0.5*(f10+f01);
        const double mlo = std::min(m1,m2), mhi = std::max(m1,m2);

        if(mhi-mlo > solid_kink*std::min(dxn(ii),dzn(kk)))
        return mlo - solid_kink_low*(mhi-mlo);

        return az*(ax*f00 + (1.0-ax)*f10) + (1.0-az)*(ax*f01 + (1.0-ax)*f11);
        }

        const double ay = wy(jj);

        const double f000 = phi(ii,jj,kk),     f100 = phi(ii+1,jj,kk),     f010 = phi(ii,jj+1,kk),     f110 = phi(ii+1,jj+1,kk);
        const double f001 = phi(ii,jj,kk+1),   f101 = phi(ii+1,jj,kk+1),   f011 = phi(ii,jj+1,kk+1),   f111 = phi(ii+1,jj+1,kk+1);
        const double m1 = 0.5*(f000+f111), m2 = 0.5*(f100+f011), m3 = 0.5*(f010+f101), m4 = 0.5*(f110+f001);
        const double mlo = std::min(std::min(m1,m2),std::min(m3,m4)), mhi = std::max(std::max(m1,m2),std::max(m3,m4));

        if(mhi-mlo > solid_kink*std::min(dxn(ii),std::min(dyn(jj),dzn(kk))))
        return mlo - solid_kink_low*(mhi-mlo);

        return      ay*(az*(ax*f000 + (1.0-ax)*f100) + (1.0-az)*(ax*f001 + (1.0-ax)*f101))
            + (1.0-ay)*(az*(ax*f010 + (1.0-ax)*f110) + (1.0-az)*(ax*f011 + (1.0-ax)*f111));
    };

    // a solid of one cell across the direction of the face (the face between the cells c and c+e,
    // one of them in the solid, the cells before and after it in the fluid)
    auto thin = [&](int ii, int jj, int kk, int di, int dj, int dk)
    {
        const bool P = phi(ii,jj,kk)<=0.0, N = phi(ii+di,jj+dj,kk+dk)<=0.0;
        if(P && !N && phi(ii-di,jj-dj,kk-dk)>0.0)
        return true;
        if(N && !P && phi(ii+2*di,jj+2*dj,kk+2*dk)>0.0)
        return true;
        return false;
    };

    // all faces between two cells of the arrays (the ghost cells as well: the cells with all faces
    // closed, solid_cell_closed)
    const int ja = jd ? j0 : 0, jb = jd ? j1 : 0;

    for(i=i0; i<i1; ++i)
    for(j=ja; j<=jb; ++j)
    for(k=k0; k<=k1; ++k)
    {
        double A;
        if(jd)
        A = solid_cut_quad(corner(i,j-1,k-1),corner(i,j,k-1),corner(i,j,k),corner(i,j-1,k));
        else
        A = solid_cut_line(corner(i,j,k-1),corner(i,j,k));

        if(thin(i,j,k,1,0,0))
        A = 0.0;

        a->apx(i,j,k) = clip(A);
    }

    if(jd)
    for(i=i0; i<=i1; ++i)
    for(j=j0; j<j1; ++j)
    for(k=k0; k<=k1; ++k)
    {
        double A = solid_cut_quad(corner(i-1,j,k-1),corner(i,j,k-1),corner(i,j,k),corner(i-1,j,k));

        if(thin(i,j,k,0,1,0))
        A = 0.0;

        a->apy(i,j,k) = clip(A);
    }

    for(i=i0; i<=i1; ++i)
    for(j=ja; j<=jb; ++j)
    for(k=k0; k<k1; ++k)
    {
        double A;
        if(jd)
        A = solid_cut_quad(corner(i-1,j-1,k),corner(i,j-1,k),corner(i,j,k),corner(i-1,j,k));
        else
        A = solid_cut_line(corner(i-1,j,k),corner(i,j,k));

        if(thin(i,j,k,0,0,1))
        A = 0.0;

        a->apz(i,j,k) = clip(A);
    }

    // the faces at the partition edges: both ranks the value of the rank of the lower cell
    pgc->start1(p,a->apx,10);
    pgc->start2(p,a->apy,11);
    pgc->start3(p,a->apz,12);
}

// D 1 2: the pressure of the cells with all faces closed (inside a solid), from the fluid.
//
// A closed cell has no row in the pressure correction; its pressure stayed that of the start.  The
// forces on a solid (P 81, the pressure a quarter of a cell off the surface, interpolated from the
// cells around it), the pressure probes next to a wall and the printed field read it.  Four layers
// of closed cells take the pressure of their neighbours towards the fluid (the mean over them),
// with the hydrostatic difference between the cell centres: constant along the normal of the wall
// for the dynamic part, exact for water at rest.  The cells deeper inside: 0.
void solid_press_extend(lexer *p, fdm *a, ghostcell *pgc, field &press, int gcval)
{
    if(p->D1!=2 || a->apx.V==nullptr || !a->ap_done)
    return;

    int i,j,k;
    const int marge = increment::marge;
    const double unset = 1.0e30;
    const bool jd = p->j_dir==1;

    std::vector<int> closed;
    LOOP
    if(solid_cell_closed(p,a,i,j,k))
    {
    closed.push_back(IJK);
    press(i,j,k) = unset;
    }

    // (level 0: all ranks together, the halo updates below are collective)
    if((p->amr_patch==0 ? pgc->globalisum((int)closed.size()) : (int)closed.size())==0)
    return;

    pgc->start4(p,press,gcval);

    std::vector<double> val;

    for(int layer=0; layer<press_layers; ++layer)
    {
        val.assign(closed.size(),unset);

        int q=0;
        LOOP
        if(solid_cell_closed(p,a,i,j,k))
        {
            if(press(i,j,k)>0.5*unset)
            {
            double sum = 0.0;
            int num = 0;

            auto take = [&](int ii, int jj, int kk, double gx, double gy, double gz)
            {
                if(p->flag4[(ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax + kk-p->kmin]<=0)
                return;
                const double pn = press(ii,jj,kk);
                if(pn>0.5*unset)
                return;
                // a patch of the mesh refinement: the cells around it are filled from the other grids,
                // a closed one of them holds no pressure of the flow
                if(p->amr_patch>0 && (ii<0 || ii>=p->knox || jj<0 || jj>=p->knoy || kk<0 || kk>=p->knoz)
                   && solid_cell_closed(p,a,ii,jj,kk))
                return;
                // hydrostatic from the neighbour to the cell (gx, gy, gz: the way times the gravity)
                sum += pn + a->ro(ii,jj,kk)*(gx+gy+gz);
                ++num;
            };

            take(i-1,j,k, p->W20*p->DXP[IM1],0.0,0.0);
            take(i+1,j,k,-p->W20*p->DXP[IP],0.0,0.0);
            if(jd)
            {
            take(i,j-1,k,0.0, p->W21*p->DYP[JM1],0.0);
            take(i,j+1,k,0.0,-p->W21*p->DYP[JP],0.0);
            }
            take(i,j,k-1,0.0,0.0, p->W22*p->DZP[KM1]);
            take(i,j,k+1,0.0,0.0,-p->W22*p->DZP[KP]);

            if(num>0)
            val[q] = sum/double(num);
            }
        ++q;
        }

        q=0;
        LOOP
        if(solid_cell_closed(p,a,i,j,k))
        {
            if(val[q]<0.5*unset)
            press(i,j,k) = val[q];
        ++q;
        }

        // the last layer: the cells left are deeper inside
        if(layer==press_layers-1)
        LOOP
        if(press(i,j,k)>0.5*unset)
        press(i,j,k) = 0.0;

        pgc->start4(p,press,gcval);
    }
}
