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

double sflow_timestep_SI::wave_cap(lexer *p, double dx, double h)
{
    const double g = fabs(p->W22);
    h = MAX(h,1.0e-6);

    // hydrostatic: the gravity wave Courant number
    if(p->A220!=2 && p->A220!=3)
    {
    const double cr = p->A228<0.0 ? -p->A228 : 10.0;
    return cr*dx/sqrt(g*h);
    }

    // non-hydrostatic: N steps per period of a wave of 10 cells, its phase speed from the dispersion of the model
    const double nt = p->A228<0.0 ? -p->A228 : 20.0;
    const double alpha = (p->A220==3) ? MAX(p->A224,1.0) : 1.0;
    const double B = (alpha-1.0)/3.0;
    const double lambda = 10.0*dx;
    const double kh = 2.0*PI/lambda*h;
    const double c = sqrt(g*h*(1.0 + B*kh*kh)/(1.0 + alpha/3.0*kh*kh));

    return lambda/(c*nt);
}

double sflow_timestep_SI::advective(lexer *p, fdm2D* b, ghostcell* pgc, double &hmax)
{
    // 1/max(|u|/dx + |v|/dy), the explicit diffusion limit (A 212 1) scaled as in
    // sflow_etimestep, and the largest water depth
    double sigma,cmin=1.0e20;
    capw = 1.0e20;

    // the non-hydrostatic pressure (A 220 2/3): the vertical momentum WH is advected with AB2 (upwind fluxes,
    // stable up to the Courant number 0.5), so the advective limit is halved; in waves |u| is small and the
    // cap A 228 decides, the limit acts in fast thin flows (run-down, bores)
    // A 210 11 (inside A210 = 1): 2.5, the flow Courant number at most 0.8 N 47 (0.4 with N 47 0.5): the one-step
    // BCG predictor misses bores from about 0.5 on (Stoker dam break, theta 1/2: max. rel. L1 0.15 % at u dt/dx 0.4,
    // 12 % at 0.5, 6 - 11 % at 1.0)
    const double nhfac = (p->A226==1 && p->A210==1) ? 2.5 : ((p->A220>=2) ? 2.0 : 1.0);

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

    sigma *= nhfac;

    if(sigma>1.0e-20)
    cmin = MIN(cmin,1.0/sigma);

    hmax = MAX(hmax,b->WL(i,j));

    // the wave limit without A 228 > 0 (the cell width: the smaller one; curvilinear grid from 1/dx^2 + 1/dy^2)
    if(p->A228<=0.0)
    {
    double dxc;

    if(p->geo_curv==1)
    dxc = sqrt(2.0/MAX(sflow_curv::inv2(p,i,j),1.0e-20));

    else
    dxc = (p->j_dir==1) ? MIN(p->DXN[IP],p->DYN[JP]) : p->DXN[IP];

    capw = MIN(capw,wave_cap(p,dxc,b->WL(i,j)));
    }
    }

    if(p->A212==1)
    SLICELOOP4
    WETDRY
    {
    double visc = p->W2 + b->eddyv(i,j);

    if(p->A246==2 && b->breaking(i,j)==1)
    visc += p->A250;

    if(p->A246==4)
    visc += b->vb(i,j);

    p->viscmax = MAX(p->viscmax,visc);

    const double lam = (p->geo_curv==1) ? sflow_curv::inv2(p,i,j) : 1.0/(p->DXN[IP]*p->DXN[IP]) + p->y_dir/(p->DYN[JP]*p->DYN[JP]);

    cmin = MIN(cmin, 0.5/(visc*lam + 1.0e-20)/(2.0*cfl));
    }

    // X 10 3: the moving pressure travels at most about one cell per step (the predictor takes it explicitly); X 10 2:
    // the body of the direct forcing
    if(p->X10==2 || p->X10==3)
    {
    const double us = fabs(p->ufbi)/p->DXM + p->y_dir*fabs(p->vfbi)/p->DXM;

    if(us>1.0e-20)
    cmin = MIN(cmin,1.0/us);
    }

    cmin = pgc->globalmin(cmin);
    hmax = pgc->globalmax(hmax);

    if(p->A228<=0.0)
    capw = pgc->globalmin(capw);
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

    // cap: A 228 > 0, otherwise the wave limit of the grid (wave_cap)
    const double cap = p->A228>0.0 ? p->A228 : capw;

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
    {
    if(p->A220==2 || p->A220==3)
    cout<<"SFLOW A 210 12/13: time step from the grid, "<<(p->A228<0.0 ? -p->A228 : 20.0)
        <<" steps per period of a wave of 10 cells (A 228 > 0: the time step in s, A 228 < 0: the steps per period)"<<endl;
    else
    cout<<"SFLOW A 210 12/13: time step from the grid, gravity wave Courant number "<<(p->A228<0.0 ? -p->A228 : 10.0)
        <<" (A 228 > 0: the time step in s, A 228 < 0: the Courant number)"<<endl;
    }

	p->dt=pgc->timesync(p->dt);
	p->dt_old=p->dt;

	b->maxF=0.0;
	b->maxG=0.0;
}
