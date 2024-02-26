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

#include"pjm_comp.h"
#include"lexer.h"
#include"fdm.h" 
#include"ghostcell.h"
#include"poisson.h"
#include"solver.h"
#include"momentum.h"
#include"ioflow.h"
#include"fluid_update_fsf_comp.h"
#include"heat.h"
#include"concentration.h"
#include"density_f.h"
#include"density_df.h"
#include"density_comp.h"
#include"density_conc.h"
#include"density_heat.h"
#include"density_vof.h"
#include"density_rheo.h"
 
pjm_comp::pjm_comp(lexer* p, fdm *a, ghostcell *pgc, heat *&pheat, concentration *&pconc) : ro_n(p)
{
    if((p->F80==0) && p->H10==0 && p->W30==0  && p->F300==0 && p->W90==0 && p->X10==0)
	pd = new density_f(p);
    
    if((p->F80==0) && p->H10==0 && p->W30==0  && p->F300==0 && p->W90==0 && p->X10==1)  
	pd = new density_df(p);
    
	if(p->F80==0 && p->H10==0 && p->W30==1  && p->F300==0 && p->W90==0)
	pd = new density_comp(p);
	
	if(p->F80==0 && p->H10>0 && p->F300==0 && p->W90==0)
	pd = new density_heat(p,pheat);
	
	if(p->F80==0 && p->C10>0 && p->F300==0 && p->W90==0)
	pd = new density_conc(p,pconc);
    
    if(p->F80>0 && p->H10==0 && p->W30==0  && p->F300==0 && p->W90==0)
	pd = new density_vof(p);
    
    if(p->F30>0 && p->H10==0 && p->W30==0  && p->F300==0 && p->W90>0)
    pd = new density_rheo(p);
    
    if(p->F300>=1)
    pd = new density_rheo(p);
    
    
    gcval_press=40;  
	
	gcval_u=7;
	gcval_v=8;
	gcval_w=9;
    
    pupdate = new fluid_update_fsf_comp(p,a,pgc);
}

pjm_comp::~pjm_comp()
{
}

void pjm_comp::start(fdm* a,lexer*p, poisson* ppois,solver* psolv, ghostcell* pgc, ioflow *pflow, field& uvel, field& vvel, field& wvel, double alpha)
{
    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<".";
			
	vel_setup(p,a,pgc,uvel,vvel,wvel,alpha);
    density_ini(p,a,pgc);	
    rhs(p,a,pgc,uvel,vvel,wvel,alpha);
	
    ppois->start(p,a,a->press);
	
        starttime=pgc->timer();

	psolv->start(p,a,pgc,a->press,a->rhsvec,5);
	
        endtime=pgc->timer();

	pgc->start4(p,a->press,gcval_press);
	
	ucorr(p,a,uvel,alpha);
	vcorr(p,a,vvel,alpha);
	wcorr(p,a,wvel,alpha);
    
    pupdate->start(p,a,pgc);
    ptimesave(p,a,pgc);

    p->poissoniter=p->solveriter;

	p->poissontime=endtime-starttime;

	if(p->mpirank==0 && (p->count%p->P12==0))
	{
	cout<<"piter: "<<p->solveriter<<"  ptime: "<<setprecision(3)<<p->poissontime;
	
	if(p->N10>10)
	cout<<"   preconiter: "<<p->preconiter;
	cout<<endl;
    }
}

void pjm_comp::ucorr(lexer* p, fdm* a, field& uvel,double alpha)
{	
	ULOOP
	uvel(i,j,k) -= alpha*p->dt*((a->press(i+1,j,k)-a->press(i,j,k))
	/(p->DXP[IP]*pd->roface(p,a,1,0,0)));
}

void pjm_comp::vcorr(lexer* p, fdm* a, field& vvel,double alpha)
{	
	VLOOP
	vvel(i,j,k) -= alpha*p->dt*((a->press(i,j+1,k)-a->press(i,j,k))
	/(p->DYP[JP]*pd->roface(p,a,0,1,0)));
}

void pjm_comp::wcorr(lexer* p, fdm* a, field& wvel,double alpha)
{	
	WLOOP
	wvel(i,j,k) -= alpha*p->dt*((a->press(i,j,k+1)-a->press(i,j,k))
	/(p->DZP[KP]*pd->roface(p,a,0,0,1)));
}

void pjm_comp::rhs(lexer *p, fdm* a, ghostcell *pgc, field &u, field &v, field &w,double alpha)
{
    double H;
    double epsi=-0.6*p->DXM;
    
    count=0;
	
    pip=p->Y50;
    LOOP
    {
        if(a->phi(i,j,k)>epsi)
		H=0.0;

		if(a->phi(i,j,k)<epsi)
		H=1.0;

		//if(fabs(a->phi(i,j,k))<=epsi)
		//H=0.5*(1.0 - a->phi(i,j,k)/epsi + (1.0/PI)*sin((PI*(-a->phi(i,j,k)))/epsi));
        
    a->rhsvec.V[count] = -((u(i,j,k)-u(i-1,j,k))
						  +(v(i,j,k)-v(i,j-1,k))
						  +(w(i,j,k)-w(i,j,k-1)) 
                         ) /(alpha*p->dt*p->DXM)
						 

                         
                         - (H*(a->ro(i,j,k)-ro_n(i,j,k)))/(a->ro(i,j,k)*alpha*p->dt)
                         /*
                         - (H*(a->u(i,j,k)*(ro_n(i+1,j,k)-ro_n(i-1,j,k))/ro_n(i,j,k)))
                         - (H*(a->v(i,j,k)*(ro_n(i,j+1,k)-ro_n(i,j-1,k))/ro_n(i,j,k)))
                         - (H*(a->w(i,j,k)*(ro_n(i,j,k+1)-ro_n(i,j,k-1))/ro_n(i,j,k))))/(alpha*p->dt)*/;
    ++count;
    }
    pip=0;
}

void pjm_comp::vel_setup(lexer *p, fdm* a, ghostcell *pgc, field &u, field &v, field &w,double alpha)
{
	pgc->start1(p,u,gcval_u);
	pgc->start2(p,v,gcval_v);
	pgc->start3(p,w,gcval_w);
}

void pjm_comp::density_ini(lexer*p,fdm* a, ghostcell *pgc)
{
    if(p->count<=1)
    {
    LOOP
    ro_n(i,j,k)=a->ro(i,j,k);
    
    pgc->start4(p,ro_n,1);
    }
}

void pjm_comp::upgrad(lexer*p,fdm* a, slice &eta, slice &eta_n)
{
}

void pjm_comp::vpgrad(lexer*p,fdm* a, slice &eta, slice &eta_n)
{
}

void pjm_comp::wpgrad(lexer*p,fdm* a, slice &eta, slice &eta_n)
{
}

void pjm_comp::ptimesave(lexer *p, fdm *a, ghostcell *pgc)
{
	LOOP
    ro_n(i,j,k)=a->ro(i,j,k);
    
    pgc->start4(p,ro_n,1);
}

void pjm_comp::ini(lexer*p,fdm* a, ghostcell *pgc)
{
}




