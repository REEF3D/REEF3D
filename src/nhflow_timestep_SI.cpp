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

#include"nhflow_timestep_SI.h"
#include<iomanip>
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"nhflow_curv.h"
#include"nhflow_momentum_SI.h"
#include"sflow_curv.h"

nhflow_timestep_SI::nhflow_timestep_SI(lexer* p) : nhflow_timestep(p), cfl(p->N47)
{
}

nhflow_timestep_SI::~nhflow_timestep_SI()
{
}

double nhflow_timestep_SI::limit(lexer *p, fdm_nhf *d, ghostcell *pgc, double &dxm, double &hmax)
{
    const double g = fabs(p->W22);
    double cmin=1.0e20;
    hmax=0.0;
    const bool hdi = nhflow_momentum_SI::hdiff_implicit(p) && (p->A512>0 || p->A560>0 || p->A547_c>0.0 || p->A547_nu0>0.0);

    p->umax=p->vmax=p->wmax=p->viscmax=p->omegamax=0.0;

    LOOP
    WETDRY
    {
	p->umax=MAX(p->umax,fabs(d->U[IJK]));
	p->vmax=MAX(p->vmax,fabs(d->V[IJK]));
	p->wmax=MAX(p->wmax,fabs(d->W[IJK]));
    }

    FLOOP
	p->omegamax=MAX(p->omegamax,fabs(d->omegaF[FIJK]));

    // advective limit per cell (horizontal; the vertical advection is implicit), largest water depth
    LOOP
    WETDRY
    if(d->WL(i,j)>p->A544)
    {
    // the speed of the first internal mode adds to the flow (A 620: the baroclinic term is explicit)
    const double ci = p->A620>0 ? nhflow_momentum_SI::internal_speed(p,d,i,j) : 0.0;

    // the Courant number of the flow and the internal speed (curvilinear grid: the speed over the widths of the cell
    // across its faces); A 549 1 (semi-Lagrangian advection): the flow up to A549_c
    cmin = MIN(cmin,nhflow_momentum_SI::adv_cmin(p,d->U[IJK],d->V[IJK],ci,i,j,cfl));

    hmax = MAX(hmax,d->WL(i,j));

        // the horizontal diffusion of the predictor is explicit (the vertical one implicit):
        // dt <= 0.25/(nu_h (1/dx^2 + 1/dy^2)); A 512 2 on one grid: implicit, no limit (viscmax for the log only)
        if(hdi)
        p->viscmax = MAX(p->viscmax, d->VISC[IJK] + (p->A561>0 ? d->EVH[IJK] : d->EV[IJK]));

        else
        if(p->A512>0 || p->A560>0 || p->A547_c>0.0 || p->A547_nu0>0.0)
        {
        const double visc = d->VISC[IJK] + (p->A561>0 ? d->EVH[IJK] : d->EV[IJK]);
        const double lam = p->geo_curv==1 ? sflow_curv::inv2(p,i,j) : 1.0/(p->DXN[IP]*p->DXN[IP]) + p->y_dir/(p->DYN[JP]*p->DYN[JP]);

        p->viscmax = MAX(p->viscmax, visc);
        cmin = MIN(cmin, 0.25/(visc*lam + 1.0e-20)/(2.0*cfl));
        }
    }

    // multi-block grid (G 50 1): the cells of the blocks
    dxm = p->DXM;

    if(phook!=nullptr)
    {
    double dxb=1.0e20, dzb=1.0e20;
    phook->dt_si_cmin(cfl,cmin,hmax);
    phook->dt_cell_size(0,dxb,dzb);
    dxm = MIN(dxm,dxb);

    p->umax=MAX(p->umax,phook->dt_local_max(1));
    p->vmax=MAX(p->vmax,phook->dt_local_max(2));
    p->wmax=MAX(p->wmax,phook->dt_local_max(3));
    p->omegamax=MAX(p->omegamax,phook->dt_local_max(4));
    }

    dxm = pgc->globalmin(dxm);

    // 6DOF (X 10 1, 2, 4; nhflow_momentum_SI_6dof.cpp): the body moves at most N 47 of a cell per step, horizontally of
    // the smallest cell size, vertically of the thinnest layer of the wet columns (the body is cast once per step)
    if(p->X10>0)
    {
        double dzm = 1.0e20;

        LOOP
        WETDRY
        if(d->WL(i,j)>p->A544)
        dzm = MIN(dzm, p->DZN[KP]*d->WL(i,j));

        dzm = pgc->globalmin(dzm);

        const double uh = MAX(p->ufbmax,p->vfbmax);
        double dtb = 1.0e20;

        if(uh>1.0e-10)
        dtb = MIN(dtb, cfl*dxm/uh);

        if(p->wfbmax>1.0e-10)
        dtb = MIN(dtb, cfl*dzm/p->wfbmax);

        cmin = MIN(cmin, dtb/(2.0*cfl));
    }

	p->umax=pgc->globalmax(p->umax);
	p->vmax=pgc->globalmax(p->vmax);
	p->wmax=pgc->globalmax(p->wmax);
	p->omegamax=pgc->globalmax(p->omegamax);
    p->viscmax=pgc->globalmax(p->viscmax);
    cmin = pgc->globalmin(cmin);
    hmax = pgc->globalmax(hmax);

    // cap: A 536, without it the explicit gravity wave limit
    const double cap = p->A536>0.0 ? p->A536 : cfl*dxm/sqrt(g*MAX(hmax,p->A544));

    return MIN(cfl*2.0*cmin, cap);
}

void nhflow_timestep_SI::start(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    const double g = fabs(p->W22);
    double dxm, hmax;

	p->dt_old=p->dt;

    // the time step grows by at most 20 % per step
    p->dt = MIN(limit(p,d,pgc,dxm,hmax), 1.2*p->dt_old);

	p->dt=pgc->timesync(p->dt);

    if(p->N48==0)
    p->dt=p->N49;

    if(p->mpirank==0 && (p->count%p->P12==0))
    {
	cout<<"umax: "<<setprecision(5)<<p->umax<<" \t utime: "<<p->utime<<endl;
	cout<<"vmax: "<<setprecision(5)<<p->vmax<<" \t vtime: "<<p->vtime<<endl;
	cout<<"wmax: "<<setprecision(5)<<p->wmax<<" \t wtime: "<<p->wtime<<endl;
    cout<<"omegamax: "<<setprecision(5)<<p->omegamax<<endl;
    cout<<"fsftime: "<<p->fsftime<<endl;
    cout<<"SI Courant number gravity waves: "<<setprecision(3)<<p->dt*sqrt(g*MAX(hmax,p->A544))/dxm<<endl;
    }

    p->recontime=0.0;
    p->fsftime=0.0;
    p->dftime=0.0;

    d->maxF=0.0;
	d->maxG=0.0;
	d->maxH=0.0;
}

void nhflow_timestep_SI::ini(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    // the explicit initial time step (gravity waves), from which the time step grows (start)
    nhflow_timestep::ini(p,d,pgc);

    if(p->A536<=0.0 && p->mpirank==0)
    cout<<"NHFLOW A 510 11/12/13: no A 536 (maximum time step) given, the time step stays at the explicit gravity wave limit"<<endl;

	p->dt=pgc->timesync(p->dt);
	p->dt_old=p->dt;
}

void nhflow_timestep_SI::ini_end(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    // the first step at the limit of the semi-implicit scheme, which has no gravity wave limit, instead of the explicit
    // initial time step (ini: the flow speed + 10 m/s plus the gravity wave speed), from which the time step needed
    // 20 to 30 steps at 20 % per step to reach A 536. All initial fields are set here (initial and potential flow,
    // inflow, wave generation zones): the advective limit sees their velocities and, as a margin, twice the inflow
    // and outflow velocities and the surface orbital velocity of the generated waves (the relaxation zones impose
    // their velocities during the first step). The turbulence model (A 560), a fixed time step (N 48 0) and 6DOF
    // (X 10) keep the explicit initial time step.
    // 6DOF (X 10): the explicit initial time step as well, the body velocity is known after the first body step
    if(p->A560!=0 || p->N48==0 || p->X10>0)
    return;

    double dxm, hmax;
    const double dtl = limit(p,d,pgc,dxm,hmax);

    double us = MAX(fabs(p->Ui),fabs(p->Uo));

    if(p->wA>0.0 && p->ww>0.0 && p->wk>0.0 && p->wd>0.0)
    us = MAX(us, p->wA*p->ww/tanh(MIN(p->wk*p->wd,20.0)));

    us = pgc->globalmax(us);

    p->dt = us>1.0e-10 ? MIN(dtl, cfl*dxm/(2.0*us)) : dtl;

    // jumps of the initial free surface (a dam break): the flow that the surface slope accelerates during the first
    // step stays within the Courant number, g (deta/dx) dt^2/dx <= N47, deta the jump between two wet cells or the
    // water depth of a wet cell next to a dry one; after the first step the advective limit sees that flow
    const double g = fabs(p->W22);
    double dtj = 1.0e20;

    auto jump = [&](int wn, double en, double dxf)
    {
        const double de = wn==1 ? fabs(d->eta(i,j) - en) : d->WL(i,j);

        if(de>1.0e-12)
        dtj = MIN(dtj, dxf*sqrt(cfl/(g*de)));
    };

    SLICELOOP4
    if(p->wet[IJ]==1)
    {
        if(p->flagslice4[Ip1J]>0)
        jump(p->wet[Ip1J],d->eta(i+1,j),p->DXP[IP]);

        if(p->flagslice4[Im1J]>0 && p->wet[Im1J]==0)
        jump(0,0.0,p->DXP[IM1]);

        if(p->j_dir==1 && p->flagslice4[IJp1]>0)
        jump(p->wet[IJp1],d->eta(i,j+1),p->DYP[JP]);

        if(p->j_dir==1 && p->flagslice4[IJm1]>0 && p->wet[IJm1]==0)
        jump(0,0.0,p->DYP[JM1]);
    }

    p->dt = MIN(p->dt, pgc->globalmin(dtj));

	p->dt=pgc->timesync(p->dt);
	p->dt_old=p->dt;

    if(p->mpirank==0)
    cout<<"NHFLOW semi-implicit: first time step "<<setprecision(4)<<p->dt<<" s"<<endl;
}
