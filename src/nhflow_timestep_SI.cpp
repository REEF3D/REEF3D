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
#include"sflow_curv.h"

nhflow_timestep_SI::nhflow_timestep_SI(lexer* p) : nhflow_timestep(p), cfl(p->N47)
{
}

nhflow_timestep_SI::~nhflow_timestep_SI()
{
}

void nhflow_timestep_SI::start(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    const double g = fabs(p->W22);
    double hmax=0.0, sigma, cmin=1.0e20;

    p->umax=p->vmax=p->wmax=p->viscmax=p->omegamax=0.0;
	p->dt_old=p->dt;

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
    sigma = fabs(d->U[IJK])/p->DXN[IP];

    if(p->j_dir==1)
    sigma += fabs(d->V[IJK])/p->DYN[JP];

    // curvilinear grid: the speed over the widths of the cell across its faces
    if(p->geo_curv==1)
    {
    const double sp = sqrt(d->U[IJK]*d->U[IJK] + d->V[IJK]*d->V[IJK]);
    sigma = sp/nhflow_curv::width1(p,i,j) + sp/nhflow_curv::width2(p,i,j);
    }

    if(sigma>1.0e-20)
    cmin = MIN(cmin,1.0/sigma);

    hmax = MAX(hmax,d->WL(i,j));

        // the horizontal diffusion of the predictor is explicit (the vertical one implicit):
        // dt <= 0.25/(nu_h (1/dx^2 + 1/dy^2))
        if(p->A512>0 || p->A560>0)
        {
        const double visc = d->VISC[IJK] + (p->A561>0 ? d->EVH[IJK] : d->EV[IJK]);
        const double lam = p->geo_curv==1 ? sflow_curv::inv2(p,i,j) : 1.0/(p->DXN[IP]*p->DXN[IP]) + p->y_dir/(p->DYN[JP]*p->DYN[JP]);

        p->viscmax = MAX(p->viscmax, visc);
        cmin = MIN(cmin, 0.25/(visc*lam + 1.0e-20)/(2.0*cfl));
        }
    }

    // multi-block grid (G 50 1): the cells of the blocks
    double dxm = p->DXM;

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

	p->umax=pgc->globalmax(p->umax);
	p->vmax=pgc->globalmax(p->vmax);
	p->wmax=pgc->globalmax(p->wmax);
	p->omegamax=pgc->globalmax(p->omegamax);
    p->viscmax=pgc->globalmax(p->viscmax);
    cmin = pgc->globalmin(cmin);
    hmax = pgc->globalmax(hmax);

    // cap: A 536, without it the explicit gravity wave limit
    const double cap = p->A536>0.0 ? p->A536 : cfl*dxm/sqrt(g*MAX(hmax,p->A544));

    p->dt = MIN(cfl*2.0*cmin, cap);

    // the time step grows by at most 20 % per step: the turbulence model (A 560) and the open
    // boundaries start from the explicit initial time step
    p->dt = MIN(p->dt, 1.2*p->dt_old);

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
    cout<<"NHFLOW A 534 1: no A 536 (maximum time step) given, the time step stays at the explicit gravity wave limit"<<endl;

	p->dt=pgc->timesync(p->dt);
	p->dt_old=p->dt;
}
