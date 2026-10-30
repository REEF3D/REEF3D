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

#include"sflow_timestep_SI.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"sflow_curv.h"

sflow_timestep_SI::sflow_timestep_SI(lexer *p, fdm2D *b)
{
    wd_criterion=p->A244;
    cfl = p->N47;
}

sflow_timestep_SI::~sflow_timestep_SI()
{
}

double sflow_timestep_SI::advective(lexer *p, fdm2D* b, ghostcell* pgc, double &hmax)
{
    // 1/max(|u|/dx + |v|/dy), the explicit diffusion limit (A 212 1) scaled as in
    // sflow_etimestep, and the largest water depth
    double sigma,cmin=1.0e20;

    hmax = 0.0;

    SLICELOOP4
    WETDRY
    {
    // curvilinear grid: the normal velocities over the face distances (sflow_curv::sigma without the wave speed)
    if(p->geo_curv==1)
    sigma = sflow_curv::sigma(p,b,i,j,0.0);

    else
    {
    sigma = fabs(b->U(i,j))/p->DXN[IP];

    if(p->j_dir==1)
    sigma += fabs(b->V(i,j))/p->DYN[JP];
    }

    if(sigma>1.0e-20)
    cmin = MIN(cmin,1.0/sigma);

    hmax = MAX(hmax,b->WL(i,j));
    }

    if(p->A212==1)
    SLICELOOP4
    WETDRY
    {
    double visc = p->W2 + b->eddyv(i,j);

    if(p->A246==2 && b->breaking(i,j)==1)
    visc += p->A250;

    p->viscmax = MAX(p->viscmax,visc);

    const double lam = (p->geo_curv==1) ? sflow_curv::inv2(p,i,j) : 1.0/(p->DXN[IP]*p->DXN[IP]) + p->y_dir/(p->DYN[JP]*p->DYN[JP]);

    cmin = MIN(cmin, 0.5/(visc*lam + 1.0e-20)/(2.0*cfl));
    }

    cmin = pgc->globalmin(cmin);
    hmax = pgc->globalmax(hmax);
    p->viscmax = pgc->globalmax(p->viscmax);

    return cmin;
}

void sflow_timestep_SI::start(lexer *p, fdm2D* b, ghostcell* pgc)
{
    const double g = fabs(p->W22);
    double cmin,hmax,dxmin;

	p->umax=p->vmax=p->viscmax=0.0;
	p->dt_old=p->dt;

	SLICELOOP4
    WETDRY
    {
	p->umax=MAX(p->umax,fabs(b->U(i,j)));
	p->vmax=MAX(p->vmax,fabs(b->V(i,j)));
    }

	p->umax=pgc->globalmax(p->umax);
	p->vmax=pgc->globalmax(p->vmax);

    cmin = advective(p,b,pgc,hmax);

    dxmin = p->DXM;

    // cap: A 228, without it the explicit gravity wave limit
    const double cap = p->A228>0.0 ? p->A228 : cfl*2.0*dxmin/sqrt(g*MAX(hmax,wd_criterion));

    p->dt = MIN(cfl*2.0*cmin, cap);

	p->dt=pgc->timesync(p->dt);

    if(p->mpirank==0 && (p->count%p->P12==0))
    {
	cout<<"umax: "<<setprecision(3)<<p->umax<<" \t utime: "<<p->utime<<endl;
	cout<<"vmax: "<<setprecision(3)<<p->vmax<<" \t vtime: "<<p->vtime<<endl;
	cout<<"fsftime: "<<p->lsmtime<<endl;
    cout<<"SI Courant number gravity waves: "<<setprecision(3)<<p->dt*sqrt(g*MAX(hmax,wd_criterion))/dxmin<<endl;
    }

	b->maxF=0.0;
	b->maxG=0.0;
}

void sflow_timestep_SI::ini(lexer *p, fdm2D* b, ghostcell* pgc)
{
    const double g = fabs(p->W22);

    // the initial velocities are not known yet: the cap A 228, otherwise the explicit
    // initial time step (gravity waves at the initial depth)
    double h = MAX(MAX(p->F60,p->wd),1.0);

    p->dt = cfl*p->DXM/(MAX(p->W10,0.0) + sqrt(g*h));

    if(p->A228>0.0)
    p->dt = p->A228;

    if(p->A228<=0.0 && p->mpirank==0)
    cout<<"SFLOW A 226 1: no A 228 (maximum time step) given, the time step stays at the explicit gravity wave limit"<<endl;

	p->dt=pgc->timesync(p->dt);
	p->dt_old=p->dt;

	b->maxF=0.0;
	b->maxG=0.0;
}
