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

#include"fem_coupling.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"slice.h"
#include<mpi.h>
#include<cmath>
#include<algorithm>
#include<map>
#include<iostream>

// FEM coupling with REEF3D::NHFLOW (sigma grid). The algorithm of start_cfd:
// direct forcing of the structural surface velocity, debris drag, parcel and
// probed pressure loads, solid time step in the final stage. Fluid data:
//   velocities U,V,W at the cell centres, the column heights bed + ZP WL
//   pressure: non-hydrostatic P on the sigma levels (bed + ZN WL), low-pass
//             filtered over a few steps, plus the hydrostatic pressure
//             W1 g (bed + WL - z) below the local free surface
//   water: wet columns below the free surface
// The kernel is the Roma kernel of CFD: horizontally over the cell centres,
// vertically over the cell centres of the wet column with the width
// max(layer thickness, smallest horizontal cell) (the Lagrangian points are
// spaced by the horizontal cells), renormalised over the column at the bed and
// the free surface. The domain is decomposed horizontally only.
// The cells inside deformable structures are solid for NHFLOW (p->DF -1, no
// fluxes through their faces). NHFLOW keeps every water column from the bed to
// the free surface: a structure that stands dry closes its columns over the full
// height, water cannot flow over it (overtopping needs CFD; a warning is given).

typedef fem_solid::Vec3 Vec3;

double fem_coupling::nhf_column(lexer *p, double *f, int ic, int jc, double z, bool faces) const
{
    // linear interpolation in one sigma column at the height z: cell centres,
    // or the sigma levels (faces); the heights are rebuilt from bed + sigma WL
    // (the ZSP / ZSN halo is not reliable at subdomain faces)
    int i = ic, j = (p->j_dir==1) ? jc : 0, k = 0;
    const double wl = (*WLn)(i,j), zb = nhf->bed(i,j);
    const int nk = faces ? p->knoz+1 : p->knoz;

    auto zl = [&](int kk)->double {return zb + (faces ? p->ZN[kk+marge] : p->ZP[kk+marge])*wl;};
    auto val = [&](int kk)->double {k = kk; return faces ? f[FIJK] : f[IJK];};

    if(nk<2 || wl<1.0e-12 || z<=zl(0))
    return val(0);

    for(int kk=0; kk<nk-1; ++kk)
    {
        const double za = zl(kk), zc = zl(kk+1);
        if(z<=zc)
        {
            const double w = (zc>za) ? (z-za)/(zc-za) : 0.0;
            return (1.0-w)*val(kk) + w*val(kk+1);
        }
    }
    return val(nk-1);
}

double fem_coupling::nhf_ipol(lexer *p, double *f, double x, double y, double z, bool faces) const
{
    // bilinear between the four columns around (x,y)
    const int i0 = p->posf_i(x);
    double wa = (p->XP[i0+1+marge]-x)/p->DXP[i0+marge];
    wa = std::max(0.0,std::min(1.0,wa));
    if(p->j_dir==0)
    return wa*nhf_column(p,f,i0,0,z,faces) + (1.0-wa)*nhf_column(p,f,i0+1,0,z,faces);

    const int j0 = p->posf_j(y);
    double wb = (p->YP[j0+1+marge]-y)/p->DYP[j0+marge];
    wb = std::max(0.0,std::min(1.0,wb));
    return wa*wb*nhf_column(p,f,i0,j0,z,faces) + (1.0-wa)*wb*nhf_column(p,f,i0+1,j0,z,faces)
         + wa*(1.0-wb)*nhf_column(p,f,i0,j0+1,z,faces) + (1.0-wa)*(1.0-wb)*nhf_column(p,f,i0+1,j0+1,z,faces);
}

double fem_coupling::nhf_bedlevel(lexer *p, double x, double y) const
{
    return p->ccslipol4(nhf->bed,x,y);
}

double fem_coupling::nhf_surface(lexer *p, double x, double y) const
{
    return p->ccslipol4(nhf->bed,x,y) + p->ccslipol4(*WLn,x,y);
}

bool fem_coupling::nhf_wet_column(lexer *p, double x, double y) const
{
    const int i = std::max(0,std::min(p->posc_i(x),p->knox-1));
    const int j = (p->j_dir==1) ? std::max(0,std::min(p->posc_j(y),p->knoy-1)) : 0;
    return p->wet[IJ]!=0 && (*WLn)(i,j)>1.0e-10;
}

double fem_coupling::nhf_level(lexer *p, const Vec3& x) const
{
    // signed distance to the bed and the immersed solids of the grid
    double lv = x(2) - nhf_bedlevel(p,x(0),x(1));
    if(nhf_solid)
    lv = std::min(lv, nhf_ipol(p,nhf->SOLID,x(0),x(1),x(2),false));
    return lv;
}

void fem_coupling::nhf_filter_pressure()
{
    // The non-hydrostatic pressure next to the forced surface jumps from step to
    // step where a body fills most of a shallow water column (the forcing of the
    // column changes the divergence the projection removes); taken as it is, it
    // drives light debris unstable. The probes take it low-pass filtered over
    // a few fluid steps; the hydrostatic part is not filtered.
    const int np = (int)pts.size();
    const double c = 1.0/3.0;
    if((int)pnh_bar.size()!=np)
    {
        pnh_bar.assign(np,0.0);
        pnh_bar_b.assign(np,0.0);
        for(int q=0; q<np; ++q)
        {
            const double *b = &buf[BP*q];
            if(b[5]>0.5) pnh_bar[q] = b[18]/b[5];
            if(b[16]>0.5) pnh_bar_b[q] = b[19]/b[16];
        }
    }
    for(int q=0; q<np; ++q)
    {
        double *b = &buf[BP*q];
        if(b[5]>0.5)
        {
            const double pn = b[18]/b[5];
            pnh_bar[q] += c*(pn - pnh_bar[q]);
            b[4] += (pnh_bar[q] - pn)*b[5];
        }
        if(b[16]>0.5)
        {
            const double pn = b[19]/b[16];
            pnh_bar_b[q] += c*(pn - pnh_bar_b[q]);
            b[15] += (pnh_bar_b[q] - pn)*b[16];
        }
    }
}

bool fem_coupling::nhf_in_solid(lexer *p, const Vec3& x) const
{
    // in the bed, a solid of the grid or a cell marked solid (p->DF < 0: the
    // structures of this coupling): no pressure probe there
    if(nhf_level(p,x)<0.0)
    return true;
    const int i = std::max(0,std::min(p->posc_i(x(0)),p->knox-1));
    const int j = (p->j_dir==1) ? std::max(0,std::min(p->posc_j(x(1)),p->knoy-1)) : 0;
    const double wl = (*WLn)(i,j), zb = nhf->bed(i,j);
    int k = 0;
    if(wl>=1.0e-10)
    while(k<p->knoz-1 && zb+p->ZN[k+1+marge]*wl<=x(2))
    ++k;
    return p->DF[IJK]<0;
}

double fem_coupling::nhf_hv(lexer *p, int ic, int jc, double z) const
{
    // vertical kernel width: the layer thickness at z, at least the smallest
    // horizontal cell (the spacing of the Lagrangian points)
    const int i = std::max(0,std::min(ic,p->knox-1));
    const int j = (p->j_dir==1) ? std::max(0,std::min(jc,p->knoy-1)) : 0;
    const double wl = (*WLn)(i,j);
    if(wl<1.0e-10)
    return dxmin;
    const double sg = (z - nhf->bed(i,j))/wl;
    int kc = 0;
    while(kc<p->knoz-1 && p->ZN[kc+1+marge]<=sg)
    ++kc;
    return std::max(dxmin, p->DZN[kc+marge]*wl);
}

double fem_coupling::nhf_kernel_ipol(lexer *p, double *f, const Vec3& x) const
{
    // kernel interpolation of a cell-centred field (owner rank, the stencil
    // reaches 2 ghost cells), wet columns only
    const int ii = p->posc_i(x(0));
    const int jj = (p->j_dir==1) ? p->posc_j(x(1)) : 0;
    const double hv = nhf_hv(p,ii,jj,x(2));

    double s = 0.0, ws = 0.0;
    for(int ci=0; ci<5; ++ci)
    {
        const int i = ii+ci-2;
        const double Dx = kernel((p->XP[i+marge]-x(0))/p->DXN[ii+marge]);
        if(Dx==0.0) continue;
        for(int cj=0; cj<5; ++cj)
        {
            if(p->j_dir==0 && cj!=2) continue;
            const int j = jj+cj-2;
            const double Dy = (p->j_dir==1) ? kernel((p->YP[j+marge]-x(1))/p->DYN[jj+marge]) : 1.0;
            if(Dy==0.0 || p->wet[IJ]==0) continue;
            const double wl = (*WLn)(i,j), zb = nhf->bed(i,j);
            if(wl<1.0e-10) continue;
            for(int k=0; k<p->knoz; ++k)
            {
                const double Dz = kernel((zb + p->ZP[k+marge]*wl - x(2))/hv);
                if(Dz==0.0) continue;
                const double D = Dx*Dy*Dz;
                s += D*f[IJK];
                ws += D;
            }
        }
    }
    return ws>0.0 ? s/ws : 0.0;
}

void fem_coupling::nhf_spread(lexer *p, double *UH, double *VH, double *WH, const Vec3& xp, const Vec3& du, double A, const Vec3* n)
{
    // velocity increment du at the point, spread onto the own cells (ghost cells
    // are exchanged by nhflow_forcing): with a normal n, A is the point area and
    // the volume A times the cell size along n; without, A is the volume.
    // Vertically renormalised over the wet column (no loss at bed and surface).
    if(xp(0)<p->originx-2.0*p->DXN[marge] || xp(0)>p->endx+2.0*p->DXN[p->knox-1+marge])
    return;
    if(p->j_dir==1 && (xp(1)<p->originy-2.0*p->DYN[marge] || xp(1)>p->endy+2.0*p->DYN[p->knoy-1+marge]))
    return;

    const int ii = p->posc_i(xp(0));
    const int jj = (p->j_dir==1) ? p->posc_j(xp(1)) : 0;
    const double hv = nhf_hv(p,ii,jj,xp(2));

    const int is = std::max(ii-2,0), ie = std::min(ii+2,p->knox-1);
    const int js = (p->j_dir==1) ? std::max(jj-2,0) : 0;
    const int je = (p->j_dir==1) ? std::min(jj+2,p->knoy-1) : 0;
    if(is>ie || js>je)
    return;

    std::vector<double> Dz(p->knoz);

    for(int i=is; i<=ie; ++i)
    {
        const double dx = p->DXN[i+marge];
        const double Dx = kernel((p->XP[i+marge]-xp(0))/dx);
        if(Dx==0.0) continue;

        for(int j=js; j<=je; ++j)
        {
            const double dy = p->DYN[j+marge];
            const double Dy = (p->j_dir==1) ? kernel((p->YP[j+marge]-xp(1))/dy) : 1.0;
            if(Dy==0.0 || p->wet[IJ]==0) continue;
            const double wl = (*WLn)(i,j), zb = nhf->bed(i,j);
            if(wl<1.0e-10) continue;

            double sum = 0.0;
            for(int k=0; k<p->knoz; ++k)
            {
                Dz[k] = kernel((zb + p->ZP[k+marge]*wl - xp(2))/hv);
                sum += Dz[k]*p->DZN[k+marge]*wl;
            }
            if(sum<1.0e-20) continue;

            double dV = A;
            if(n)
            dV = A*(std::fabs((*n)(0))*dx + std::fabs((*n)(1))*(p->j_dir==1 ? dy : 0.0) + std::fabs((*n)(2))*hv);
            // 2D: forces are per unit width, the fluid cell has width dy
            if(p->j_dir==0)
            dV *= dy;

            const double fac = dV*(Dx/dx)*(Dy/dy)/sum;

            for(int k=0; k<p->knoz; ++k)
            {
                if(Dz[k]==0.0) continue;
                const double w = fac*Dz[k];
                nhf->U[IJK] += w*du(0);
                UH[IJK]     += w*du(0)*wl;
                if(p->j_dir==1)
                {
                nhf->V[IJK] += w*du(1);
                VH[IJK]     += w*du(1)*wl;
                }
                nhf->W[IJK] += w*du(2);
                WH[IJK]     += w*du(2)*wl;
            }
        }
    }
}

void fem_coupling::nhf_probe_pressure(lexer *p, const Vec3& xp, const Vec3& n, double *b)
{
    // pressure probe outside the surface point (owner rank of the probe):
    // b[4] pressure, b[5] count, b[6] water at the probe; the offset is
    // shortened when the probe would leave the domain horizontally
    double off = fs.coupling().pressure_offset;
    auto in_domain = [&](const Vec3& q)
    {
        const double eps = 1.0e-9;
        return q(0)>p->global_xmin+eps && q(0)<p->global_xmax-eps
            && (p->j_dir==0 || (q(1)>p->global_ymin+eps && q(1)<p->global_ymax-eps));
    };
    const int ic = p->posc_i(xp(0)), jc = (p->j_dir==1) ? p->posc_j(xp(1)) : 0;
    const int ib = std::max(0,std::min(ic,p->knox-1)), jb = std::max(0,std::min(jc,p->knoy-1));
    const Vec3 h(p->DXN[ib+marge], p->j_dir==1 ? p->DYN[jb+marge] : 0.0, nhf_hv(p,ic,jc,xp(2)));

    for(double f : {1.0, 2.0/3.0, 1.0/3.0})
    {
        off = f*fs.coupling().pressure_offset;
        if(in_domain(xp + 1.01*off*n.cwiseProduct(h)))
        break;
    }
    Vec3 pr = xp + off*n.cwiseProduct(h);
    if(p->j_dir==0)
    pr(1) = p->YP[marge];

    auto owns = [&](const Vec3& q){return q(0)>=p->originx && q(0)<p->endx && (p->j_dir==0 || (q(1)>=p->originy && q(1)<p->endy));};
    if(!owns(pr))
    return;

    // a probe in the bed or a solid is moved closer to the surface, then out horizontally
    auto inside = [&](const Vec3& q){return nhf_in_solid(p,q);};
    for(double f : {2.0/3.0, 1.0/3.0})
    if(inside(pr))
    {
        Vec3 q = xp + f*off*n.cwiseProduct(h);
        if(p->j_dir==0) q(1) = pr(1);
        if(!inside(q))
        pr = q;
    }
    if(inside(pr))
    {
        Vec3 nh(n(0),n(1),0.0);
        if(nh.norm()>0.3)
        {
            nh.normalize();
            pr = xp + off*nh.cwiseProduct(h);
            if(p->j_dir==0) pr(1) = p->YP[marge];
        }
        if(!owns(pr))
        return;
    }
    if(inside(pr))
    return;

    b[5] = 1.0;
    if(!nhf_wet_column(p,pr(0),pr(1)))
    return;     // dry: no pressure, no water

    // non-hydrostatic pressure at the probe, hydrostatic pressure of the local
    // free surface at the height of the surface point
    const double eta = nhf_surface(p,pr(0),pr(1));
    const double pnh = nhf_ipol(p,nhf->P,pr(0),pr(1),pr(2),true);
    b[4] = pnh + p->W1*std::fabs(p->W22)*std::max(0.0,eta-xp(2));
    b[18] = pnh;
    b[6] = pr(2)<=eta ? 1.0 : 0.0;
}

void fem_coupling::nhf_probe_beside(lexer *p, const Vec3& xp, Vec3 q, double *b)
{
    // pressure beside a rigid body at the height of xp: b[0] pressure, b[1] count, b[2] water
    if(p->j_dir==0)
    q(1) = p->YP[marge];
    if(q(0)<p->originx || q(0)>=p->endx)
    return;
    if(p->j_dir==1 && (q(1)<p->originy || q(1)>=p->endy))
    return;
    q(2) = std::max(q(2), nhf_bedlevel(p,q(0),q(1)));
    for(int lift=0; lift<2 && nhf_in_solid(p,q); ++lift)
    q(2) += 0.5*dxmin;
    if(nhf_in_solid(p,q))
    return;

    b[1] = 1.0;
    if(!nhf_wet_column(p,q(0),q(1)))
    return;
    const double eta = nhf_surface(p,q(0),q(1));
    const double pnh = nhf_ipol(p,nhf->P,q(0),q(1),q(2),true);
    b[0] = pnh + p->W1*std::fabs(p->W22)*std::max(0.0,eta-xp(2));
    b[4] = pnh;     // [19] of the point
    b[2] = xp(2)<=eta ? 1.0 : 0.0;
}

void fem_coupling::nhf_sample_bed(lexer *p, ghostcell *pgc)
{
    // signed distance to the bed and the solids of the grid and its gradient at
    // the nodes of free bodies and debris near the bed; owner rank, one reduction
    const int nn = fs.nnode();
    std::vector<double> b(5*size_t(nn),0.0);
    const double hb = 2.0*std::max(fs.hmin(),dxmin);

    for(int i=0; i<nn; ++i)
    {
        if(!fs.free_node(i))
        continue;
        Vec3 x = fs.pos(i);
        if(p->j_dir==0)
        x(1) = p->YP[marge];
        if(x(0)<p->originx || x(0)>=p->endx)
        continue;
        if(p->j_dir==1 && (x(1)<p->originy || x(1)>=p->endy))
        continue;

        const double phi = nhf_level(p,x);
        if(phi>hb)
        continue;

        const double d = 0.5*dxmin;
        Vec3 g;
        g(0) = (nhf_level(p,x+Vec3(d,0.0,0.0)) - nhf_level(p,x-Vec3(d,0.0,0.0)))/(2.0*d);
        g(1) = (p->j_dir==1) ? (nhf_level(p,x+Vec3(0.0,d,0.0)) - nhf_level(p,x-Vec3(0.0,d,0.0)))/(2.0*d) : 0.0;
        g(2) = (nhf_level(p,x+Vec3(0.0,0.0,d)) - nhf_level(p,x-Vec3(0.0,0.0,d)))/(2.0*d);
        double* q = &b[5*size_t(i)];
        q[0] = phi; q[1] = g(0); q[2] = g(1); q[3] = g(2); q[4] = 1.0;
    }

    if(!b.empty())
    MPI_Allreduce(MPI_IN_PLACE,b.data(),(int)b.size(),MPI_DOUBLE,MPI_SUM,pgc->mpi_comm);

    fs.clear_bed_samples();
    for(int i=0; i<nn; ++i)
    {
        const double* q = &b[5*size_t(i)];
        if(q[4]<0.5)
        continue;
        Vec3 n(q[1]/q[4],q[2]/q[4],q[3]/q[4]);
        const double gn = n.norm();
        if(gn<1.0e-6)
        continue;
        fs.set_bed_sample(i,q[0]/q[4],n/gn);
    }
}

void fem_coupling::nhf_mark_solid(lexer *p, ghostcell *pgc, bool dfreset)
{
    // The kernel forcing alone does not stop the flow through a wall in NHFLOW:
    // the Riemann fluxes between the columns carry water into and through it.
    // As for the solids of the grid, the cells inside the intact elements of the
    // deformable bodies get p->DF = -1 (wall states at their faces, and no
    // continuity flux, d->solid_flux): cells whose centre lies in the bounding box
    // of an element, and the cell of every element centre (structures thinner
    // than a cell). A dry column under a structure has all its cells at the bed:
    // it is closed over the full height. Rigid bodies (debris) are not marked.
    // The marks of the previous stage are taken back unless nhflow_forcing has
    // reset p->DF.
    if(!dfreset)
    for(size_t q=0; q<df_cells.size(); ++q)
    if(p->DF[df_cells[q]]==-1)
    p->DF[df_cells[q]] = df_save[q];
    df_cells.clear();
    df_save.clear();

    // dry columns covered by the structure from the bed: top of the structure,
    // for the overtopping warning
    std::map<int,double> dry_top;

    auto mark = [&](int i, int j, int k)
    {
        const int c = IJK;
        if(p->DF[c]==-1)
        return;
        df_cells.push_back(c);
        df_save.push_back(p->DF[c]);
        p->DF[c] = -1;
    };

    for(int e=0; e<fs.nelem(); ++e)
    {
        const fem_solid::element& el = fs.elem(e);
        if(!el.alive || el.rigid || fs.rigid_of_node(el.n[0])>=0)
        continue;

        Vec3 lo = fs.pos(el.n[0]), hi = lo, c = Vec3::Zero();
        for(int a=0; a<8; ++a)
        {
            lo = lo.cwiseMin(fs.pos(el.n[a]));
            hi = hi.cwiseMax(fs.pos(el.n[a]));
            c += 0.125*fs.pos(el.n[a]);
        }
        if(hi(0)<p->originx-p->DXN[marge] || lo(0)>p->endx+p->DXN[p->knox-1+marge])
        continue;
        if(p->j_dir==1 && (hi(1)<p->originy-p->DYN[marge] || lo(1)>p->endy+p->DYN[p->knoy-1+marge]))
        continue;

        // cell centres inside the bounding box
        const int is = std::max(0,p->posc_i(lo(0))-1), ie = std::min(p->knox-1,p->posc_i(hi(0))+1);
        const int js = (p->j_dir==1) ? std::max(0,p->posc_j(lo(1))-1) : 0;
        const int je = (p->j_dir==1) ? std::min(p->knoy-1,p->posc_j(hi(1))+1) : 0;
        for(int i=is; i<=ie; ++i)
        {
            if(p->XP[i+marge]<lo(0) || p->XP[i+marge]>hi(0)) continue;
            for(int j=js; j<=je; ++j)
            {
                if(p->j_dir==1 && (p->YP[j+marge]<lo(1) || p->YP[j+marge]>hi(1))) continue;
                const double wl = (*WLn)(i,j), zb = nhf->bed(i,j);
                for(int k=0; k<p->knoz; ++k)
                {
                    const double zc = zb + p->ZP[k+marge]*wl;
                    if(zc>=lo(2) && zc<=hi(2))
                    mark(i,j,k);
                }
                if(p->wet[IJ]==0 && hi(2)>zb)
                {
                    auto it = dry_top.find(IJ);
                    if(it==dry_top.end()) dry_top[IJ] = hi(2);
                    else it->second = std::max(it->second,hi(2));
                }
            }
        }

        // the cell of the element centre
        if(c(0)>=p->originx && c(0)<p->endx && (p->j_dir==0 || (c(1)>=p->originy && c(1)<p->endy)))
        {
            const int i = std::max(0,std::min(p->posc_i(c(0)),p->knox-1));
            const int j = (p->j_dir==1) ? std::max(0,std::min(p->posc_j(c(1)),p->knoy-1)) : 0;
            const double wl = (*WLn)(i,j), zb = nhf->bed(i,j);
            if(c(2)>=zb && (wl<1.0e-10 || c(2)<=zb+wl))
            {
                int k = 0;
                if(wl>=1.0e-10)
                while(k<p->knoz-1 && zb+p->ZN[k+1+marge]*wl<=c(2))
                ++k;
                mark(i,j,k);
            }
        }
    }

    pgc->startintV(p,p->DF,1);

    // NHFLOW keeps the water columns from the bed to the free surface: a dry
    // column inside a structure stays closed over its full height, the water
    // cannot flow over an emerged structure (as for the solids of NHFLOW)
    if(!warned_overtop)
    {
        int over = 0;
        for(const auto& c : dry_top)
        {
            const int i = c.first/p->jmax + p->imin, j = c.first%p->jmax + p->jmin;
            const int di[4] = {1,-1,0,0}, dj[4] = {0,0,1,-1};
            for(int q=0; q<(p->j_dir==1 ? 4 : 2); ++q)
            {
                const int in = i+di[q], jn = j+dj[q];
                if(p->wet[(in-p->imin)*p->jmax + (jn-p->jmin)]!=0 && nhf->bed(in,jn)+(*WLn)(in,jn) > c.second + 0.5*dxmin)
                over = 1;
            }
        }
        if(pgc->globalimax(over)>0)
        {
            warned_overtop = true;
            if(p->mpirank==0)
            std::cout<<"FEM WARNING: NHFLOW: the water beside the structure is higher than the structure (t = "<<p->simtime
                     <<" s): NHFLOW does not model overtopping, the structure holds the water back over the full depth "
                       "(water flows over a structure only where the structure stood in water from the start); use REEF3D::CFD for overtopping"<<std::endl;
        }
    }
}

void fem_coupling::start_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, double alpha,
                                double *UH, double *VH, double *WH, slice &WL, bool solid, bool dfreset, bool finalize)
{
    starttime = pgc->timer();
    nhf = d;
    WLn = &WL;
    nhf_solid = solid;
    d->solid_flux = 1;

    if(!initialised)
    first_call(p,pgc);

    if(surf_version!=fs.surface_version())
    ini_points(p,pgc);

    // the structure as a solid for the fluxes of the next stage
    nhf_mark_solid(p,pgc,dfreset);

    // halo of the stage velocity for the kernel interpolation
    pgc->start4V(p,d->U,10);
    pgc->start4V(p,d->V,11);
    pgc->start4V(p,d->W,12);

    const int np = (int)pts.size();
    const std::vector<int>& deb = fs.debris();
    const int nd = (int)deb.size();
    const int nn = fs.nnode();
    const bool parcels = (fs.coupling().loads!=1);
    const bool probes = (fs.coupling().loads!=2) || fs.n_rigid()>0;

    const size_t nbuf = BP*size_t(np+nd) + (finalize ? 2*size_t(nn) : 0);
    if(buf.size()<nbuf)
    buf.resize(nbuf);
    std::fill(buf.begin(),buf.begin()+nbuf,0.0);

    auto owns = [&](const Vec3& x)->bool
    {
        if(x(0)<p->originx || x(0)>=p->endx)
        return false;
        if(p->j_dir==1 && (x(1)<p->originy || x(1)>=p->endy))
        return false;
        return true;
    };
    // water at x: wet column, between the bed and the free surface
    auto in_water = [&](const Vec3& x, double& eta)->bool
    {
        eta = nhf_surface(p,x(0),x(1));
        return nhf_wet_column(p,x(0),x(1)) && x(2)<=eta && x(2)>=nhf_bedlevel(p,x(0),x(1));
    };

    // ------------------------------------------------------------------
    // 1. sample (owner rank), buffer layout as in start_cfd:
    //    [0-3] velocity, count  [4-6] pressure, count, water at probe
    //    [7] cell size normal to the surface  [8-10] momentum  [11-13] density
    //    [14] distance below the free surface  [15-17] probe beside rigid bodies
    //    points above the free surface are not sampled (no forcing, no parcel)
    // ------------------------------------------------------------------
    Vec3 xp, vp, n;
    double A, eta;

    if(finalize && probes && fs.n_rigid()>0)
    rigid_boxes();

    for(int q=0; q<np; ++q)
    {
        point_state(q,xp,vp,n,A);
        double *b = &buf[BP*q];
        Vec3 x2 = xp;
        if(p->j_dir==0)
        x2(1) = p->YP[marge];

        if(owns(x2) && in_water(x2,eta))
        {
            b[0] = nhf_kernel_ipol(p,d->U,x2);
            b[1] = (p->j_dir==1) ? nhf_kernel_ipol(p,d->V,x2) : 0.0;
            b[2] = nhf_kernel_ipol(p,d->W,x2);
            b[3] = 1.0;
            b[14] = eta - xp(2);

            if(finalize && parcels)
            {
                const int ii = std::max(0,std::min(p->posc_i(xp(0)),p->knox-1));
                const int jj = (p->j_dir==1) ? std::max(0,std::min(p->posc_j(xp(1)),p->knoy-1)) : 0;
                b[7] = std::fabs(n(0))*p->DXN[ii+marge] + std::fabs(n(1))*(p->j_dir==1 ? p->DYN[jj+marge] : 0.0) + std::fabs(n(2))*nhf_hv(p,ii,jj,xp(2));
                b[8] = p->W1*b[0];  b[11] = p->W1;
                if(p->j_dir==1)
                {
                    b[9] = p->W1*b[1];  b[12] = p->W1;
                }
                b[10] = p->W1*b[2];  b[13] = p->W1;
            }
        }

        if(finalize && probes)
        {
            nhf_probe_pressure(p,xp,n,b);
            Vec3 fb;
            if(probe_fallback(q,xp,fb))
            nhf_probe_beside(p,xp,fb,b+15);
        }
    }

    // debris particles: velocity, count, water, density
    for(int dd=0; dd<nd; ++dd)
    {
        Vec3 x = fs.pos(deb[dd]);
        if(p->j_dir==0)
        x(1) = p->YP[marge];
        double *b = &buf[BP*(np+dd)];

        if(owns(x))
        {
            const bool w = in_water(x,eta);
            b[0] = w ? nhf_ipol(p,d->U,x(0),x(1),x(2),false) : 0.0;
            b[1] = (w && p->j_dir==1) ? nhf_ipol(p,d->V,x(0),x(1),x(2),false) : 0.0;
            b[2] = w ? nhf_ipol(p,d->W,x(0),x(1),x(2),false) : 0.0;
            b[3] = 1.0;
            b[4] = w ? 1.0 : 0.0;
            b[5] = p->W1;
        }
    }

    // fluid density at the nodes (enclosed fluid of the immersed boundary)
    if(finalize && parcels)
    for(int i=0; i<nn; ++i)
    {
        Vec3 x = fs.pos(i);
        if(p->j_dir==0)
        x(1) = p->YP[marge];
        if(owns(x))
        {
            double *b = &buf[BP*size_t(np+nd)+2*size_t(i)];
            b[0] = in_water(x,eta) ? p->W1 : 0.0;
            b[1] = 1.0;
        }
    }

    if(nbuf>0)
    MPI_Allreduce(MPI_IN_PLACE,buf.data(),(int)nbuf,MPI_DOUBLE,MPI_SUM,pgc->mpi_comm);

    // ------------------------------------------------------------------
    // 2. direct forcing of the structural surface velocity: the velocity
    //    increment (u_s - u) is spread onto U,V,W and UH,VH,WH
    // ------------------------------------------------------------------
    const double hy = fs.lattice_h(1);

    if(fs.coupling().forcing)
    for(int q=0; q<np; ++q)
    {
        const double *b = &buf[BP*q];
        if(b[3]<0.5)
        continue;

        point_state(q,xp,vp,n,A);

        if(on_bed(q,n,p->global_zmin))
        continue;

        Vec3 uf(b[0]/b[3],b[1]/b[3],b[2]/b[3]);
        if(p->j_dir==0)
        uf(1) = vp(1);

        const Vec3 du = vp - uf;
        const double Aeff = (p->j_dir==0) ? A/hy : A;
        Vec3 x2 = xp;
        if(p->j_dir==0)
        x2(1) = p->YP[marge];
        nhf_spread(p,UH,VH,WH,x2,du,Aeff,&n);
    }

    // ------------------------------------------------------------------
    // 3. debris drag (quadratic), reaction onto the fluid
    // ------------------------------------------------------------------
    const double cd = fs.coupling().debris_cd;

    for(int dd=0; dd<nd; ++dd)
    {
        const double *b = &buf[BP*(np+dd)];
        fdeb[dd].setZero();

        if(b[3]<0.5 || b[4]/b[3]<0.5)
        continue;

        const int i = deb[dd];
        Vec3 ur = Vec3(b[0],b[1],b[2])/b[3] - fs.vel(i);
        if(p->j_dir==0)
        ur(1) = 0.0;

        const double rho = b[5]/b[3];
        if(!ur.allFinite() || !std::isfinite(rho) || rho<=0.0)
        continue;
        const double Ad = std::pow(fs.node_volume(i),2.0/3.0);
        fdeb[dd] = 0.5*rho*cd*Ad*ur.norm()*ur;

        if(fs.coupling().debris_reaction)
        {
            Vec3 fr = -fdeb[dd]/rho;
            if(p->j_dir==0)
            fr /= hy;
            Vec3 x = fs.pos(i);
            if(p->j_dir==0)
            x(1) = p->YP[marge];
            nhf_spread(p,UH,VH,WH,x,alpha*p->dt*fr,1.0,nullptr);
        }
    }

    // ghost cells of U,V,W,UH,VH,WH are updated by nhflow_forcing::forcing

    // ------------------------------------------------------------------
    // 4. loads and solid time step (final stage)
    // ------------------------------------------------------------------
    if(finalize)
    {
        if(fs.bed_contact())
        nhf_sample_bed(p,pgc);
        if(probes)
        nhf_filter_pressure();
        finish_step(p,pgc,alpha);
    }
}
