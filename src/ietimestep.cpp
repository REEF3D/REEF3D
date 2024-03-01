/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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
#include"ietimestep.h"
#include<iomanip>
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"turbulence.h"

ietimestep::ietimestep(lexer* p):epsi(1.0e-19),c0_orig(p->N47)
{
}

ietimestep::~ietimestep()
{
}

void ietimestep::start(fdm *a, lexer *p, ghostcell *pgc, turbulence *pturb)
{
    p->umax=p->vmax=p->wmax=p->viscmax=irsm=jrsm=krsm=0.0;
    p->epsmax=p->kinmax=p->pressmax=0.0;
	p->pressmin=1.0e9;

	p->umax=p->vmax=p->wmax=p->viscmax=0.0;
    
    p->umax=MAX(p->W11_u,p->umax);
    p->umax=MAX(p->W12_u,p->umax);
    p->umax=MAX(p->W13_u,p->umax);
    p->umax=MAX(p->W14_u,p->umax);
    p->umax=MAX(p->W15_u,p->umax);
    p->umax=MAX(p->W16_u,p->umax);
    
    p->vmax=MAX(p->W11_v,p->vmax);
    p->vmax=MAX(p->W12_v,p->vmax);
    p->vmax=MAX(p->W13_v,p->vmax);
    p->vmax=MAX(p->W14_v,p->vmax);
    p->vmax=MAX(p->W15_v,p->vmax);
    p->vmax=MAX(p->W16_v,p->vmax);
    
    p->wmax=MAX(p->W11_w,p->wmax);
    p->wmax=MAX(p->W12_w,p->wmax);
    p->wmax=MAX(p->W13_w,p->wmax);
    p->wmax=MAX(p->W14_w,p->wmax);
    p->wmax=MAX(p->W15_w,p->wmax);
    p->wmax=MAX(p->W16_w,p->wmax);
    
	sqd=1.0/(p->DXM*p->DXM);
	
// maximum velocities

	ULOOP
	p->umax=MAX(p->umax,fabs(a->u(i,j,k)));

	p->umax=pgc->globalmax(p->umax);
    
    


	VLOOP
	p->vmax=MAX(p->vmax,fabs(a->v(i,j,k)));

	p->vmax=pgc->globalmax(p->vmax);


	WLOOP
	p->wmax=MAX(p->wmax,fabs(a->w(i,j,k)));

	p->wmax=pgc->globalmax(p->wmax);
    
    velmax=max(p->umax,p->vmax,p->wmax);

    if(p->mpirank==0 && (p->count%p->P12==0))
    {
	cout<<"umax: "<<setprecision(3)<<p->umax<<" \t utime: "<<p->utime<<endl;
	cout<<"vmax: "<<setprecision(3)<<p->vmax<<" \t vtime: "<<p->vtime<<endl;
	cout<<"wmax: "<<setprecision(3)<<p->wmax<<" \t wtime: "<<p->wtime<<endl;
    }
	
	p->umax=MAX(p->umax,p->ufbmax);
	p->vmax=MAX(p->vmax,p->vfbmax);
	p->wmax=MAX(p->wmax,p->wfbmax);

    velmax=max(p->umax,p->vmax,p->wmax);
    
    // rhs globalmax
    a->maxF=pgc->globalmax(a->maxF);
    a->maxG=pgc->globalmax(a->maxG);
    a->maxH=pgc->globalmax(a->maxH);
    p->fbmax=pgc->globalmax(p->fbmax);
    p->fbmax=pgc->globalmax(p->fbmax);

    // maximum viscosity
	LOOP
	p->viscmax=MAX(p->viscmax, a->visc(i,j,k)+a->eddyv(i,j,k));

	p->viscmax=pgc->globalmax(p->viscmax);

    if(p->mpirank==0 && (p->count%p->P12==0))
	cout<<"viscmax: "<<p->viscmax<<endl;
    
	//----kin
	LOOP
	p->kinmax=MAX(p->kinmax,pturb->kinval(i,j,k));

	p->kinmax=pgc->globalmax(p->kinmax);

    if(p->mpirank==0 && (p->count%p->P12==0))
	cout<<"kinmax: "<<p->kinmax<<endl;

	//---eps
    LOOP
	p->epsmax=MAX(p->epsmax,pturb->epsval(i,j,k));

	p->epsmax=pgc->globalmax(p->epsmax);

    if(p->mpirank==0 && (p->count%p->P12==0))
	cout<<"epsmax: "<<p->epsmax<<endl;


	//---press
    LOOP
    {
	p->pressmax=MAX(p->pressmax,a->press(i,j,k));
	p->pressmin=MIN(p->pressmin,a->press(i,j,k));
    }

	p->pressmax=pgc->globalmax(p->pressmax);
	p->pressmin=pgc->globalmin(p->pressmin);


// maximum reynolds stress source term
	visccrit=p->viscmax*(6.0/pow(p->DXM,2.0));
	
    
    cu=1.0e10;
    cv=1.0e10;
    cw=1.0e10;
    cb=1.0e10;
    
    if(p->N50==1)
    LOOP
    {
    dx = MIN3(p->DXN[IP],p->DYN[JP],p->DZN[KP]);

	cu = MIN(cu, 2.0/((sqrt(p->umax*p->umax + p->vmax*p->vmax + p->wmax*p->wmax))/dx
    
            + sqrt((4.0*fabs(MAX3(a->maxF,a->maxG,a->maxH)))/dx)));
    }
    
    if(p->N50==2)
    LOOP
    {
	cu = MIN(cu, 2.0/((sqrt(p->umax*p->umax))/p->DXN[IP]
    
            + sqrt((4.0*fabs(a->maxF))/p->DXN[IP])));
            
    cv = MIN(cv, 2.0/((sqrt(p->vmax*p->vmax))/p->DYN[JP]
    
            + sqrt((4.0*fabs(a->maxG))/p->DYN[JP])));
            
    cw = MIN(cw, 2.0/((sqrt(p->wmax*p->wmax))/p->DZN[KP]
    
            + sqrt((4.0*fabs(a->maxH))/p->DZN[KP])));
    }
    
    cu = MIN3(cu,cv,cw);
    
	p->dt=p->N47*cu;
	p->dt=pgc->timesync(p->dt);
    
    
    // fbdt
    LOOP
    {
    dx = MIN3(p->DXN[IP],p->DYN[JP],p->DZN[KP]);

	cb = MIN(cb, 2.0/sqrt((4.0*fabs(p->fbmax))/dx));
    }
    
    p->fbdt=p->N47*cb;
    p->fbdt=pgc->timesync(p->fbdt);

	a->maxF=0.0;
	a->maxG=0.0;
	a->maxH=0.0;
    
}

void ietimestep::ini(fdm* a, lexer* p,ghostcell* pgc)
{  
    dx = p->DXM;
    
	p->umax=p->vmax=p->wmax=p->viscmax=-1e19;
    
    p->viscmax = MAX(p->W2,p->W4);

	p->umax=MAX(p->W10,p->umax);
    
    p->umax=MAX(p->W11_u,p->umax);
    p->umax=MAX(p->W12_u,p->umax);
    p->umax=MAX(p->W13_u,p->umax);
    p->umax=MAX(p->W14_u,p->umax);
    p->umax=MAX(p->W15_u,p->umax);
    p->umax=MAX(p->W16_u,p->umax);
    
    p->vmax=MAX(p->W11_v,p->vmax);
    p->vmax=MAX(p->W12_v,p->vmax);
    p->vmax=MAX(p->W13_v,p->vmax);
    p->vmax=MAX(p->W14_v,p->vmax);
    p->vmax=MAX(p->W15_v,p->vmax);
    p->vmax=MAX(p->W16_v,p->vmax);
    
    p->wmax=MAX(p->W11_w,p->wmax);
    p->wmax=MAX(p->W12_w,p->wmax);
    p->wmax=MAX(p->W13_w,p->wmax);
    p->wmax=MAX(p->W14_w,p->wmax);
    p->wmax=MAX(p->W15_w,p->wmax);
    p->wmax=MAX(p->W16_w,p->wmax);
    
    ULOOP
	p->umax=MAX(p->umax,fabs(a->u(i,j,k)));

	p->umax=pgc->globalmax(p->umax);


	VLOOP
	p->vmax=MAX(p->vmax,fabs(a->v(i,j,k)));

	p->vmax=pgc->globalmax(p->vmax);


	WLOOP
	p->wmax=MAX(p->wmax,fabs(a->w(i,j,k)));

	p->wmax=pgc->globalmax(p->wmax);

	p->umax=MAX(p->umax,2.0*p->ufbmax);
	p->umax=MAX(p->umax,2.0*p->vfbmax);
	p->umax=MAX(p->umax,2.0*p->wfbmax);
    
    p->umax=MAX(p->umax,2.0*p->X210_u);
	p->umax=MAX(p->umax,2.0*p->X210_v);
	p->umax=MAX(p->umax,2.0*p->X210_w);
    
    p->umax=MAX(p->umax,2.0);
    
    
    cu=1.0e10;
    cv=1.0e10;
    cw=1.0e10;
    
    if(p->N50==1)
    LOOP
    {
    dx = MIN3(p->DXN[IP],p->DYN[JP],p->DZN[KP]);

	cu = MIN(cu, 2.0/((sqrt(p->umax*p->umax + p->vmax*p->vmax + p->wmax*p->wmax))/dx
    
            + sqrt((4.0*fabs(MAX3(a->maxF,a->maxG,a->maxH)))/dx)));
    }
    
    if(p->N50==2)
    LOOP
    {
	cu = MIN(cu, 2.0/((sqrt(p->umax*p->umax))/p->DXN[IP]
    
            + sqrt((4.0*fabs(a->maxF))/p->DXN[IP])));
            
    cv = MIN(cv, 2.0/((sqrt(p->vmax*p->vmax))/p->DYN[JP]
    
            + sqrt((4.0*fabs(a->maxG))/p->DYN[JP])));
            
    cw = MIN(cw, 2.0/((sqrt(p->wmax*p->wmax))/p->DZN[KP]
    
            + sqrt((4.0*fabs(a->maxH))/p->DZN[KP])));
    }
    
    cu = MIN3(cu,cv,cw);
    
	p->dt=p->N47*cu*0.25;
    p->dt = MAX(p->dt,1.0e-6);
    
	p->dt=pgc->timesync(p->dt);
	p->dt_old=p->dt;
    
    a->maxF = fabs(a->gi);
    a->maxG = fabs(a->gj);
    a->maxH = fabs(a->gk);
    
}

double ietimestep::min(double val1,double val2,double val3)
{
	double mini;

	mini=val1;

	if(mini>val2)
	mini=val2;

	if(mini>val3)
	mini=val3;

	if(mini<0.0)
	mini=0.0;

	return mini;
}

double ietimestep::min(double val1,double val2)
{
	double mini;

	mini=val1;

	if(mini>val2)
	mini=val2;

	if(mini<0.0)
	mini=0.0;

	return mini;
}

double ietimestep::max(double val1,double val2,double val3)
{
	double maxi;

	maxi=val1;

	if(maxi<val2)
	maxi=val2;

	if(maxi<val3)
	maxi=val3;

	if(maxi<0.0)
	maxi=0.0;

	return maxi;
}

double ietimestep::max(double val1,double val2)
{
	double maxi;

	maxi=val1;

	if(maxi<val2)
	maxi=val2;

	if(maxi<0.0)
	maxi=0.0;

	return maxi;
}

