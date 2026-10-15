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
#include"field.h"
#include"ioflow.h"
#include"wave_interface.h"
#include"bc_zone.h"
#include"definitions.h"
#include<cmath>
#include<iostream>

//  The relaxation zones of the wave generation (B 98 2) and of the numerical beach (B 99 1/2) on
//  the CFD patches.  Without it the refinement keeps out of the zones (reefamr3d::mark_tiles) and
//  the waves are made on level 0, whose step has to resolve them as well (G 7 1, validation/22).
//
//  A patch in a zone relaxes as level 0 does (iowave_relax.cpp, iowave_precalc_relax.cpp): the
//  velocities and the pressure after each projection, the level set after its transport (before
//  the reinitialisation), with the weights of the zones (bc_zone_set, the same zones as iowave) at
//  the positions of the patch's faces and cells, the targets of the wave theory of level 0
//  (wave_interface) evaluated there at the start of the patch's step (G 7 1; G 7 0: of the step,
//  as level 0), the interface band of the patch (3 and 0.6 of its cell size).  For wave theories
//  (not the wavemakers B 92 20-29), B 95 0, without zone sources or backgrounds (B 523, B 524);
//  otherwise the zones stay on level 0 as before.  The eddy viscosity of a RANS model relaxes in
//  the generation zone as on level 0 (iowave::turb_relax: towards 0 with the weight of the zone).

namespace
{
inline int cix_r(const lexer *q, int i, int j, int k)
{
    return (i-q->imin)*q->jmax*q->kmax + (j-q->jmin)*q->kmax + k-q->kmin;
}

// the Heaviside function of the relaxation (iowave_relax.cpp): H = G
inline double relax_H(double phival, double epsi, double psi)
{
    if(phival>=-psi)
    return 1.0;
    if(phival<-epsi)
    return 0.0;
    return 0.5*(1.0 + phival/epsi + (1.0/PI)*sin((PI*phival)/epsi));
}

// the wave frame (iowave_dist.cpp)
inline double xframe(double x, double y, double x0, double y0, double g) { return (x-x0)*cos(g) + (y-y0)*sin(g); }
inline double yframe(double x, double y, double x0, double y0, double g) { return -(x-x0)*sin(g) + (y-y0)*cos(g); }

// the wave theory at the time t of the patch (wave_interface reads it from the lexer)
struct wtime
{
    lexer *p;
    double st, wt;
    wtime(lexer *pp, double t) : p(pp), st(pp->simtime), wt(pp->wavetime) { p->simtime = t; p->wavetime = t; }
    ~wtime() { p->simtime = st; p->wavetime = wt; }
};
}

void cfd_amr::relax_setup(lexer *p, ghostcell *pgc)
{
    pwave0 = dynamic_cast<wave_interface*>(pflow0);
    if(pwave0==nullptr || !(p->B98==2 || p->B99==1 || p->B99==2))
    return;

    if(p->B95!=0.0 || (p->B92>=20 && p->B92<=29))
    {
        if(p->mpirank==0)
        cout<<"CFD AMR: the relaxation zones stay on level 0 (B 95 or a wavemaker B 92 20-29)"<<endl;
        return;
    }

    zones0 = new bc_zone_set(bc_zone_set::from_legacy(p,pgc));
    if(zones0->has_sources() || zones0->has_background() || !zones0->edges.empty())
    {
        if(p->mpirank==0)
        cout<<"CFD AMR: the relaxation zones stay on level 0 (zone sources, backgrounds or edges, B 520-524)"<<endl;
        delete zones0;
        zones0 = nullptr;
        return;
    }

    zgen = (p->B98==2) ? 1 : 0;
    zbeach = (p->B99==1 || p->B99==2) ? 1 : 0;
    zone_gen_refine = zgen==1;
    zone_beach_refine = zbeach==1;
}

// a cell of the patch (interior or the cells around it) in a zone
bool cfd_amr::patch_in_zone(int id)
{
    if(zones0==nullptr)
    return false;
    lexer *pp = CP(id)->pp;
    for(int i=-1; i<=pp->knox; ++i)
    for(int j=0; j<pp->knoy; ++j)
    {
        const double x = pp->XP[i+marge], y = pp->YP[j+marge];
        if((zgen && zones0->relax_dist(x,y)<1.0e20) || (zbeach && zones0->beach_dist(x,y)<1.0e19))
        return true;
    }
    return false;
}

// the velocities (after the projection) and the pressure of patch id in the zones
void cfd_amr::patch_relax_vel(int id, field &u, field &v, field &w)
{
    if(zones0==nullptr || !patch_in_zone(id))
    return;

    cfd_amr_patch *c = CP(id);
    lexer *pp = c->pp;
    fdm *a = c->a;
    const double t = sub ? pp->simtime : p0->simtime;
    wtime wt(p0,t);

    const double epsi = 3.0*pp->DXM, psi = 0.6*pp->DXM;
    const double gam = p0->B105_1*(PI/180.0), x0 = p0->B105_2, y0 = p0->B105_3;
    const double phimean = p0->phimean;
    const int nx = pp->knox, ny = pp->knoy, nz = pp->knoz;

    // ramp (iowave::ramp, B 95 0)
    double ramp = 1.0;
    if(p0->B101==1 && t<p0->B102*p0->wT)
    ramp = t/(p0->B102*p0->wT) - (1.0/PI)*sin(PI*(t/(p0->B102*p0->wT)));
    if(p0->B101==2 && t<p0->B102)
    ramp = t/p0->B102 - (1.0/PI)*sin(PI*(t/p0->B102));

    // the surface elevation of the wave theory at the cell centres, columns -1 .. n (generation zone)
    const int sx = nx+2, sy = ny+2;
    vector<double> eta((size_t)sx*sy,0.0);
    auto E = [&](int i, int j) -> double& { return eta[(size_t)(i+1)*sy + (j+1)]; };
    for(int i=-1; i<=nx; ++i)
    for(int j=-1; j<=ny; ++j)
    {
        const int jj = pp->j_dir==1 ? j : 0;
        const double x = pp->XP[i+marge], y = pp->YP[jj+marge];
        if(zgen && zones0->relax_dist(x,y)<1.0e20)
        E(i,j) = pwave0->wave_eta(p0,pgc0,xframe(x,y,x0,y0,gam),yframe(x,y,x0,y0,gam));
    }

    // u
    for(int i=0; i<nx; ++i)
    for(int j=0; j<ny; ++j)
    {
        const double xc = pp->XP[i+marge], yc = pp->YP[j+marge];
        const double dg = zgen ? zones0->relax_dist(xc,yc) : 1.0e20;
        const double db = zbeach ? zones0->beach_dist(xc,yc) : 1.0e20;
        if(dg>=1.0e20 && db>=1.0e19)
        continue;
        const double xu = pp->XN[i+1+marge];
        const double wg = (dg<1.0e20) ? zones0->relax_weight(xu,yc) : 1.0;
        const double wb = zones0->beach_weight(xu,yc);
        const double et = 0.5*(E(i,j)+E(i+1,j));
        for(int k=0; k<nz; ++k)
        {
            if(pp->flag1[cix_r(pp,i,j,k)]<=0)
            continue;
            const double H = relax_H(0.5*(a->phi(i,j,k)+a->phi(i+1,j,k)),epsi,psi);
            if(dg<1.0e20)
            {
                const double zl = pp->ZP[k+marge], fsf = et + phimean;
                const double z = (zl<=fsf) ? zl-phimean : et;
                const double uval = (zl<=fsf+epsi) ? pwave0->wave_u(p0,pgc0,xframe(xu,yc,x0,y0,gam),yframe(xu,yc,x0,y0,gam),z) + p0->Ui : 0.0;
                double &U = u(i,j,k);
                U = (1.0-wg)*ramp*uval*H + wg*H*U + (1.0-H)*U;
                U += (1.0-wg)*ramp*p0->W50*(1.0-H) + wg*(1.0-H)*p0->W50;
            }
            if(db<1.0e19)
            u(i,j,k) = wb*u(i,j,k);
        }
    }

    // v (3D)
    if(pp->j_dir==1)
    for(int i=0; i<nx; ++i)
    for(int j=0; j<ny; ++j)
    {
        const double xc = pp->XP[i+marge], yc = pp->YP[j+marge];
        const double dg = zgen ? zones0->relax_dist(xc,yc) : 1.0e20;
        const double db = zbeach ? zones0->beach_dist(xc,yc) : 1.0e20;
        if(dg>=1.0e20 && db>=1.0e20)
        continue;
        const double yv = pp->YN[j+1+marge];
        const double wg = (dg<1.0e20) ? zones0->relax_weight(xc,yv) : 1.0;
        const double wb = zones0->beach_weight(xc,yv);
        const double et = 0.5*(E(i,j)+E(i,j+1));
        for(int k=0; k<nz; ++k)
        {
            if(pp->flag2[cix_r(pp,i,j,k)]<=0)
            continue;
            const double H = relax_H(0.5*(a->phi(i,j,k)+a->phi(i,j+1,k)),epsi,psi);
            if(dg<1.0e20)
            {
                const double zl = pp->ZP[k+marge], fsf = et + phimean;
                const double z = (zl<=fsf) ? zl-phimean : et;
                const double vval = (zl<=fsf+epsi) ? pwave0->wave_v(p0,pgc0,xframe(xc,yv,x0,y0,gam),yframe(xc,yv,x0,y0,gam),z) : 0.0;
                double &V = v(i,j,k);
                V = (1.0-wg)*ramp*vval*H + wg*H*V + (1.0-H)*V;
            }
            if(db<1.0e20)
            v(i,j,k) = wb*v(i,j,k);
        }
    }

    // w and the pressure (cell columns)
    for(int i=0; i<nx; ++i)
    for(int j=0; j<ny; ++j)
    {
        const double xc = pp->XP[i+marge], yc = pp->YP[j+marge];
        const double dg = zgen ? zones0->relax_dist(xc,yc) : 1.0e20;
        const double db = zbeach ? zones0->beach_dist(xc,yc) : 1.0e20;
        if(dg>=1.0e20 && db>=1.0e20)
        continue;
        const double wg = (dg<1.0e20) ? zones0->relax_weight(xc,yc) : 1.0;
        const double wb = zones0->beach_weight(xc,yc);
        const double et = E(i,j);
        for(int k=0; k<nz; ++k)
        {
            const int q = cix_r(pp,i,j,k);
            if(pp->flag3[q]>0)
            {
                const double H = relax_H(0.5*(a->phi(i,j,k)+a->phi(i,j,k+1)),epsi,psi);
                if(dg<1.0e20)
                {
                    const double zl = pp->ZN[k+1+marge], fsf = et + phimean;
                    const double z = (zl<=fsf) ? zl-phimean : et;
                    const double wval = (zl<=fsf+epsi) ? pwave0->wave_w(p0,pgc0,xframe(xc,yc,x0,y0,gam),yframe(xc,yc,x0,y0,gam),z) : 0.0;
                    double &W = w(i,j,k);
                    W = (1.0-wg)*ramp*wval*H + wg*H*W + (1.0-H)*W;
                }
                if(db<1.0e20)
                w(i,j,k) = wb*w(i,j,k);
            }
            if(pp->flag4[q]>0 && db<1.0e20)
            a->press(i,j,k) = (1.0-wb)*((phimean - pp->ZP[k+marge])*a->ro(i,j,k)*fabs(p0->W22)) + wb*a->press(i,j,k);
        }
    }
}

// the level set of patch id in the zones (after its transport, before the reinitialisation)
void cfd_amr::patch_relax_phi(int id, field &f)
{
    if(zones0==nullptr || !patch_in_zone(id))
    return;

    cfd_amr_patch *c = CP(id);
    lexer *pp = c->pp;
    const double t = sub ? pp->simtime : p0->simtime;
    wtime wt(p0,t);

    const double gam = p0->B105_1*(PI/180.0), x0 = p0->B105_2, y0 = p0->B105_3;
    const double phimean = p0->phimean;

    double ramp = 1.0;
    if(p0->B101==1 && t<p0->B102*p0->wT)
    ramp = t/(p0->B102*p0->wT) - (1.0/PI)*sin(PI*(t/(p0->B102*p0->wT)));
    if(p0->B101==2 && t<p0->B102)
    ramp = t/p0->B102 - (1.0/PI)*sin(PI*(t/p0->B102));

    for(int i=0; i<pp->knox; ++i)
    for(int j=0; j<pp->knoy; ++j)
    {
        const double xc = pp->XP[i+marge], yc = pp->YP[j+marge];
        const double dg = zgen ? zones0->relax_dist(xc,yc) : 1.0e20;
        const double db = zbeach ? zones0->beach_dist(xc,yc) : 1.0e20;
        if(dg>=1.0e20 && db>=1.0e20)
        continue;
        const double wg = (dg<1.0e20) ? zones0->relax_weight(xc,yc) : 1.0;
        const double wb = zones0->beach_weight(xc,yc);
        const double et = (dg<1.0e20) ? pwave0->wave_eta(p0,pgc0,xframe(xc,yc,x0,y0,gam),yframe(xc,yc,x0,y0,gam)) : 0.0;
        for(int k=0; k<pp->knoz; ++k)
        {
            if(pp->flag4[cix_r(pp,i,j,k)]<=0)
            continue;
            const double zc = pp->ZP[k+marge];
            if(dg<1.0e20)
            f(i,j,k) = (1.0-wg)*ramp*(et + phimean - zc) + wg*f(i,j,k);
            if(db<1.0e20)
            f(i,j,k) = (1.0-wb)*(phimean - zc) + wb*f(i,j,k);
        }
    }
}

// the eddy viscosity of patch id in the generation zone (iowave::turb_relax); true if any cell
bool cfd_amr::patch_relax_turb(int id)
{
    if(zones0==nullptr || !zgen || !patch_in_zone(id))
    return false;

    cfd_amr_patch *c = CP(id);
    lexer *pp = c->pp;
    fdm *a = c->a;
    bool any = false;
    for(int i=0; i<pp->knox; ++i)
    for(int j=0; j<pp->knoy; ++j)
    {
        const double xc = pp->XP[i+marge], yc = pp->YP[j+marge];
        if(zones0->relax_dist(xc,yc)>=1.0e20)
        continue;
        const double wg = zones0->relax_weight(xc,yc);
        for(int k=0; k<pp->knoz; ++k)
        if(pp->flag4[cix_r(pp,i,j,k)]>0)
        {
            a->eddyv(i,j,k) = wg*a->eddyv(i,j,k);
            any = true;
        }
    }
    return any;
}
