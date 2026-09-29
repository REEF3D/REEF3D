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

#include"nhflow_particle_f.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"slice.h"
#include<cmath>
#include<algorithm>

// time integration: Heun (explicit trapezoidal) for the deterministic drift
//   x* = x^n + dt*u(x^n,t^n)                  (u(x^n,t^n) stored in step_begin)
//   x^n+1 = x^n + dt/2*(u(x^n,t^n) + u(x*,t^n+1))
// followed by the random walk (Euler-Maruyama, Visser 1997 drift correction for K(z)).
// In NHFLOW the waves are phase-resolved, so Stokes drift comes out of the tracking itself;
// no Stokes drift parameterisation is added.

void nhflow_particle_f::step_begin(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    seed(p,d,pgc,p->simtime+p->dt,false);

    for(auto &a : P)
    drift_velocity(p,d,a,a.u0,a.v0,a.w0);
}

void nhflow_particle_f::step_end(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    const double dt = p->dt;
    const double sqh = sqrt(2.0*std::max(p->L31,0.0)*dt);

    for(auto &a : P)
    {
        if(a.state==NHFP_STRANDED || a.state==NHFP_BED)
        continue;

        const double xo=a.x, yo=a.y, zo=a.z;
        const double u0=a.u0, v0=a.v0, w0=a.w0;

        // predictor
        a.x = xo + dt*u0;
        a.y = yo + dt*v0;
        a.z = zo + dt*w0;

        double u1,v1,w1;
        drift_velocity(p,d,a,u1,v1,w1);

        // corrector
        a.u = 0.5*(u0+u1);
        a.v = 0.5*(v0+v1);
        a.w = 0.5*(w0+w1);

        a.x = xo + dt*a.u;
        a.y = yo + dt*a.v;
        a.z = zo + dt*a.w;

        // horizontal random walk
        if(sqh>0.0)
        {
        a.x += sqh*gauss(rng);

        if(p->j_dir==1)
        a.y += sqh*gauss(rng);
        }
    }

    domain_boundary(p);

    xchange(p,pgc);

    // column dependent steps, done on the owning rank
    for(auto &a : P)
    {
        if(a.state==NHFP_BED)
        continue;

        if(a.state!=NHFP_STRANDED)
        {
        vertical_diffusion(p,d,a);
        vertical_bounds(p,d,a);
        }

        wetdry(p,d,a);
    }

    print(p,d,pgc);
}

// deterministic particle velocity: fluid + rise velocity + windage at the surface
void nhflow_particle_f::drift_velocity(lexer *p, fdm_nhf *d, nhflow_particle_data &a, double &up, double &vp, double &wp)
{
    up=vp=wp=0.0;

    if(a.state==NHFP_STRANDED || a.state==NHFP_BED)
    return;

    const bool atsurf = (a.mode>0.5 || a.state==NHFP_SURFACE);

    double zp = atsurf ? surface(p,d,a.x,a.y) : a.z;

    fluid_velocity(p,d,a.x,a.y,zp,up,vp,wp);

    if(atsurf)
    {
    up += a.cw*U10*cosw;
    vp += a.cw*U10*sinw;
    }

    if(a.mode>0.5)
    wp = 0.0; // surface-trapped: z follows eta in vertical_bounds
    else
    wp += a.ws;

    if(p->j_dir==0)
    vp=0.0;
}

// Velocity on the sigma grid: bilinear in the horizontal over the surrounding cell-centre columns,
// and linear in the vertical within each column using that column's own ZSP levels. The vertical ghost
// cells are never used. Between the top cell centre and the free surface the velocity is extrapolated
// linearly from the two top cell centres, so floating particles see the surface velocity (the top cell
// centre alone underestimates the Stokes drift by about exp(-2k*dz_top/2)). Below the lowest cell centre
// the value is held constant. Dry columns are left out of the horizontal weights.
void nhflow_particle_f::fluid_velocity(lexer *p, fdm_nhf *d, double xp, double yp, double zp, double &up, double &vp, double &wp)
{
    up=vp=wp=0.0;

    const int ic = p->posf_i(xp);
    const int jc = (p->j_dir==1) ? p->posf_j(yp) : 0;

    double wa = (p->XP[ic+1+marge]-xp)/p->DXP[ic+marge];
    wa = std::max(0.0,std::min(1.0,wa));

    double wb = 1.0;
    if(p->j_dir==1)
    {
    wb = (p->YP[jc+1+marge]-yp)/p->DYP[jc+marge];
    wb = std::max(0.0,std::min(1.0,wb));
    }

    const int nj = (p->j_dir==1) ? 2 : 1;
    double wsum=0.0;

    for(int qi=0; qi<2; ++qi)
    for(int qj=0; qj<nj; ++qj)
    {
        const int ii = ic+qi;
        const int jj = jc+qj;
        const double w = (qi==0 ? wa : 1.0-wa)*(qj==0 ? wb : 1.0-wb);

        if(w<=0.0)
        continue;

        if(p->wet[(ii-p->imin)*p->jmax + (jj-p->jmin)]<=0)
        continue;

        double uc,vc,wc;
        column_velocity(p,d,ii,jj,zp,uc,vc,wc);

        up += w*uc;
        vp += w*vc;
        wp += w*wc;
        wsum += w;
    }

    if(wsum>1.0e-12)
    {
    up/=wsum;
    vp/=wsum;
    wp/=wsum;
    }

    if(p->j_dir==0)
    vp=0.0;

    if(up!=up || vp!=vp || wp!=wp)
    up=vp=wp=0.0;
}

void nhflow_particle_f::column_velocity(lexer *p, fdm_nhf *d, int ii, int jj, double zp, double &uc, double &vc, double &wc)
{
    const int nz = p->knoz;
    const int base = (ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax - p->kmin;

    auto zc = [&](int kk){ return p->ZSP[base+kk]; };

    if(nz==1 || zp<=zc(0))
    {
    const int id0 = base;
    uc = d->U[id0];
    vc = (p->j_dir==1) ? d->V[id0] : 0.0;
    wc = d->W[id0];
    return;
    }

    int k0;
    if(zp>=zc(nz-1))
    k0 = nz-2;  // linear extrapolation towards the free surface
    else
    {
    int lo=0, hi=nz-1;
    while(hi-lo>1)
    {
        int mid=(lo+hi)/2;
        if(zc(mid)<=zp)
        lo=mid;
        else
        hi=mid;
    }
    k0=lo;
    }

    const double dz = zc(k0+1)-zc(k0);
    const double f = (dz>1.0e-12) ? (zp-zc(k0))/dz : 0.0;

    const int id0 = base+k0;
    const int id1 = base+k0+1;

    uc = d->U[id0] + f*(d->U[id1]-d->U[id0]);
    vc = (p->j_dir==1) ? d->V[id0] + f*(d->V[id1]-d->V[id0]) : 0.0;
    wc = d->W[id0] + f*(d->W[id1]-d->W[id0]);
}

// scalar on the sigma grid (eddy viscosity), same stencil as fluid_velocity, constant beyond the outer cell centres
double nhflow_particle_f::scalar_ipol(lexer *p, double *f, double xp, double yp, double zp)
{
    const int ic = p->posf_i(xp);
    const int jc = (p->j_dir==1) ? p->posf_j(yp) : 0;

    double wa = std::max(0.0,std::min(1.0,(p->XP[ic+1+marge]-xp)/p->DXP[ic+marge]));
    double wb = (p->j_dir==1) ? std::max(0.0,std::min(1.0,(p->YP[jc+1+marge]-yp)/p->DYP[jc+marge])) : 1.0;

    const int nj = (p->j_dir==1) ? 2 : 1;
    double val=0.0, wsum=0.0;

    for(int qi=0; qi<2; ++qi)
    for(int qj=0; qj<nj; ++qj)
    {
        const int ii = ic+qi;
        const int jj = jc+qj;
        const double w = (qi==0 ? wa : 1.0-wa)*(qj==0 ? wb : 1.0-wb);

        if(w<=0.0 || p->wet[(ii-p->imin)*p->jmax + (jj-p->jmin)]<=0)
        continue;

        val += w*column_scalar(p,f,ii,jj,zp);
        wsum += w;
    }

    return (wsum>1.0e-12) ? val/wsum : 0.0;
}

double nhflow_particle_f::column_scalar(lexer *p, double *f, int ii, int jj, double zp)
{
    const int nz = p->knoz;
    const int base = (ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax - p->kmin;

    if(nz==1 || zp<=p->ZSP[base])
    return f[base];

    if(zp>=p->ZSP[base+nz-1])
    return f[base+nz-1];

    int lo=0, hi=nz-1;
    while(hi-lo>1)
    {
        int mid=(lo+hi)/2;
        if(p->ZSP[base+mid]<=zp)
        lo=mid;
        else
        hi=mid;
    }

    const double dz = p->ZSP[base+lo+1]-p->ZSP[base+lo];
    const double fr = (dz>1.0e-12) ? (zp-p->ZSP[base+lo])/dz : 0.0;

    return f[base+lo] + fr*(f[base+lo+1]-f[base+lo]);
}

double nhflow_particle_f::surface(lexer *p, fdm_nhf *d, double xp, double yp)
{
    return p->ccslipol4(d->WL,xp,yp) + p->ccslipol4(d->bed,xp,yp);
}

double nhflow_particle_f::bedlevel(lexer *p, fdm_nhf *d, double xp, double yp)
{
    return p->ccslipol4(d->bed,xp,yp);
}

// K(z) = background Kv + nu_t/Sc
double nhflow_particle_f::diffusivity(lexer *p, fdm_nhf *d, double xp, double yp, double zp)
{
    double K = std::max(p->L33,0.0);

    if(p->L32==2)
    K += std::max(scalar_ipol(p,d->EV,xp,yp,zp),0.0)/std::max(p->L34,1.0e-6);

    return K;
}

void nhflow_particle_f::vertical_diffusion(lexer *p, fdm_nhf *d, nhflow_particle_data &a)
{
    if(p->L32==0 || a.mode>0.5)
    return;

    const double dt = p->dt;

    if(p->L32==1)
    {
    a.z += sqrt(2.0*std::max(p->L33,0.0)*dt)*gauss(rng);
    return;
    }

    // Visser (1997): z^n+1 = z + K'(z) dt + sqrt(2 K(z + K'dt/2) dt) N(0,1)
    double zb = bedlevel(p,d,a.x,a.y);
    double zs = surface(p,d,a.x,a.y);
    double h = std::max(zs-zb,1.0e-6);
    double dz = std::max(0.25*h/double(std::max(p->knoz,1)),1.0e-6);

    double z1 = std::max(a.z-dz,zb);
    double z2 = std::min(a.z+dz,zs);
    double dKdz = (z2>z1) ? (diffusivity(p,d,a.x,a.y,z2)-diffusivity(p,d,a.x,a.y,z1))/(z2-z1) : 0.0;

    double zm = std::min(std::max(a.z + 0.5*dKdz*dt,zb),zs);
    double K = diffusivity(p,d,a.x,a.y,zm);

    a.z += dKdz*dt + sqrt(2.0*K*dt)*gauss(rng);
}

// free surface and bed
void nhflow_particle_f::vertical_bounds(lexer *p, fdm_nhf *d, nhflow_particle_data &a)
{
    double zb = bedlevel(p,d,a.x,a.y);
    double zs = surface(p,d,a.x,a.y);

    if(a.mode>0.5)
    {
    a.z = zs;
    a.state = NHFP_SURFACE;
    return;
    }

    if(a.z>=zs)
    {
    a.z = zs;
    a.state = NHFP_SURFACE;
    return;
    }

    if(a.z<=zb)
    {
        if(p->L52==1 && a.ws<0.0)
        {
        a.z = zb;
        a.state = NHFP_BED;
        a.u=a.v=a.w=0.0;
        return;
        }

        a.z = std::min(2.0*zb-a.z,zs);
    }

    a.state = NHFP_WATER;
}

// wetting and drying: stranding in dry cells, refloating (L 51 = 0) when the cell is wet again
void nhflow_particle_f::wetdry(lexer *p, fdm_nhf *d, nhflow_particle_data &a)
{
    if(a.state==NHFP_BED)
    return;

    int ii=i, jj=j;

    i = p->posc_i(a.x);
    i = std::max(0,std::min(i,p->knox-1));

    j = 0;
    if(p->j_dir==1)
    {
    j = p->posc_j(a.y);
    j = std::max(0,std::min(j,p->knoy-1));
    }

    const bool wet = (p->wet[IJ]>0);

    i=ii;
    j=jj;

    if(!wet)
    {
    a.state = NHFP_STRANDED;
    a.z = bedlevel(p,d,a.x,a.y);
    a.u=a.v=a.w=0.0;
    return;
    }

    if(a.state==NHFP_STRANDED && p->L51==0)
    {
    a.state = NHFP_WATER;
    a.z = surface(p,d,a.x,a.y);  // refloat at the surface
    vertical_bounds(p,d,a);
    }
}

// lateral boundaries of the global domain: L 53 = 0 remove, 1 reflect
void nhflow_particle_f::domain_boundary(lexer *p)
{
    const double xs=p->global_xmin, xe=p->global_xmax;
    const double ys=p->global_ymin, ye=p->global_ymax;

    size_t q=0;
    while(q<P.size())
    {
        nhflow_particle_data &a = P[q];

        bool out = (a.x<xs || a.x>xe);

        if(p->j_dir==1)
        out = out || (a.y<ys || a.y>ye);

        if(out && p->L53==1)
        {
            if(a.x<xs) a.x = 2.0*xs-a.x;
            if(a.x>xe) a.x = 2.0*xe-a.x;
            a.x = std::min(std::max(a.x,xs),xe);

            if(p->j_dir==1)
            {
            if(a.y<ys) a.y = 2.0*ys-a.y;
            if(a.y>ye) a.y = 2.0*ye-a.y;
            a.y = std::min(std::max(a.y,ys),ye);
            }

            out=false;
        }

        if(out)
        {
            P[q] = P.back();
            P.pop_back();
            ++numout;
            continue;
        }

        ++q;
    }
}
