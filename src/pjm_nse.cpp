/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2022 Hans Bihs

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
#include"pjm_nse.h"
#include"lexer.h"
#include"fdm.h" 
#include"ghostcell.h"
#include"poisson.h"
#include"solver.h"
#include"momentum.h"
#include"ioflow.h"
#include"heat.h"
#include"concentration.h"
#include"density_f.h"#include"density_df.h"
#include"density_comp.h"
#include"density_conc.h"
#include"density_heat.h"
#include"density_vof.h"
#include"density_rheo.h"
 
pjm_nse::pjm_nse(lexer* p, fdm *a, heat *&pheat, concentration *&pconc) 
{
	if((p->F80==0||p->A10==5) && p->H10==0 && p->W30==0  && p->F300==0 && p->W90==0 && p->X10==0)	pd = new density_f(p);        if((p->F80==0||p->A10==5) && p->H10==0 && p->W30==0  && p->F300==0 && p->W90==0 && p->X10==1)  	pd = new density_df(p);
    
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
    
    if(p->F300==1)
    pd = new density_rheo(p);
    

    gcval_press=40;  
	
	gcval_u=7;
	gcval_v=8;
	gcval_w=9;
}

pjm_nse::~pjm_nse()
{
}

void pjm_nse::start(fdm* a,lexer*p, poisson* ppois,solver* psolv, ghostcell* pgc, ioflow *pflow, field& uvel, field& vvel, field& wvel, double alpha)
{
    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<".";
			
	vel_setup(p,a,pgc,uvel,vvel,wvel,alpha);	
    rhs(p,a,pgc,uvel,vvel,wvel,alpha);
	
    
    pgc->start4(p,a->press,gcval_press);
    
    ppois->start(p,a,a->press);
	
        starttime=pgc->timer();

    psolv->start(p,a,pgc,a->press,a->rhsvec,5);
	
        endtime=pgc->timer();

	pgc->start4(p,a->press,gcval_press);
	
	ucorr(p,a,uvel,alpha);
	vcorr(p,a,vvel,alpha);
	wcorr(p,a,wvel,alpha);

    p->poissoniter=p->solveriter;

	p->poissontime=endtime-starttime;

	if(p->mpirank==0 && (p->count%p->P12==0))
	cout<<"piter: "<<p->solveriter<<"  ptime: "<<setprecision(3)<<p->poissontime<<endl;
}

void pjm_nse::ucorr(lexer* p, fdm* a, field& uvel,double alpha)
{	
    if(p->D37==1)
	ULOOP
	uvel(i,j,k) -= alpha*p->dt*CPOR1*PORVAL1*((a->press(i+1,j,k)-a->press(i,j,k))
	/(p->DXP[IP]*pd->roface(p,a,1,0,0)));
    
    if(p->D37==2)
	ULOOP
    {
    check=0;
        
        if(p->D37==2)
        {
        if(p->flag4[Ip1JK]==AIR)
        {
        teta = fabs(a->phi(i,j,k))/(fabs(a->phi(i+1,j,k))+fabs(a->phi(i,j,k))) + 0.0001*p->DXN[IP]/(fabs(a->phi(i,j,k+1))+fabs(a->phi(i,j,k)));
        check=1;
        }
        
        if(p->flag4[Im1JK]==AIR)
        {
        teta = fabs(a->phi(i,j,k))/(fabs(a->phi(i-1,j,k))+fabs(a->phi(i,j,k))) + 0.0001*p->DXN[IP]/(fabs(a->phi(i,j,k-1))+fabs(a->phi(i,j,k)));
        check=2;
        }
        
        if(p->flag4[IJp1K]==AIR)
        {
        teta = fabs(a->phi(i,j,k))/(fabs(a->phi(i,j+1,k))+fabs(a->phi(i,j,k))) + 0.0001*p->DYN[JP]/(fabs(a->phi(i,j+1,k))+fabs(a->phi(i,j,k)));
        check=1;
        }
        
        if(p->flag4[IJm1K]==AIR)
        {
        teta = fabs(a->phi(i,j,k))/(fabs(a->phi(i,j-1,k))+fabs(a->phi(i,j,k))) + 0.0001*p->DYN[JP]/(fabs(a->phi(i,j-1,k))+fabs(a->phi(i,j,k)));
        check=2;
        }
    
        if(p->flag4[IJKm1]==AIR)
        {
        teta = fabs(a->phi(i,j,k))/(fabs(a->phi(i,j,k-1))+fabs(a->phi(i,j,k))) + 0.0001*p->DZN[KP]/(fabs(a->phi(i,j,k-1))+fabs(a->phi(i,j,k)));
        check=2;
        }
        }
        
        
        if(p->flag4[IJKp1]==AIR)
        {
        teta = fabs(a->phi(i,j,k))/(fabs(a->phi(i,j,k+1))+fabs(a->phi(i,j,k))) + 0.0001*p->DZN[KP]/(fabs(a->phi(i,j,k+1))+fabs(a->phi(i,j,k)));
        check=1;
        }
        
        
    if(check==1)    
    uvel(i,j,k) -= alpha*p->dt*CPOR1*PORVAL1*(((1.0 - 1.0/teta)*a->press(i,j,k)-a->press(i,j,k))
    /(p->DXP[IP]*pd->roface(p,a,1,0,0)));
    
    if(check==2)    
    uvel(i,j,k) -= alpha*p->dt*CPOR1*PORVAL1*((a->press(i,j,k)-(1.0/teta)*a->press(i,j,k))
    /(p->DXP[IP]*pd->roface(p,a,1,0,0)));

        
    if(check==0)
    uvel(i,j,k) -= alpha*p->dt*CPOR1*PORVAL1*((a->press(i+1,j,k)-a->press(i,j,k))
	/(p->DXP[IP]*pd->roface(p,a,1,0,0)));
    }
}

void pjm_nse::vcorr(lexer* p, fdm* a, field& vvel,double alpha)
{	
	VLOOP
	vvel(i,j,k) -= alpha*p->dt*CPOR2*PORVAL2*((a->press(i,j+1,k)-a->press(i,j,k))
	/(p->DYP[JP]*pd->roface(p,a,0,1,0)));
}

void pjm_nse::wcorr(lexer* p, fdm* a, field& wvel,double alpha)
{	    
    if(p->D37==1)
	WLOOP
	wvel(i,j,k) -= alpha*p->dt*CPOR3*PORVAL3*((a->press(i,j,k+1)-a->press(i,j,k))
	/(p->DZP[KP]*pd->roface(p,a,0,0,1)));
    
    if(p->D37==2)
	WLOOP
    {
    check=0;
    
        if(p->D37==3)
        {
        if(p->flag4[Ip1JK]==AIR)
        {
        teta = fabs(a->phi(i,j,k))/(fabs(a->phi(i+1,j,k))+fabs(a->phi(i,j,k))) + 0.0001*p->DXN[IP]/(fabs(a->phi(i,j,k+1))+fabs(a->phi(i,j,k)));
        check=1;
        }
        
        if(p->flag4[Im1JK]==AIR)
        {
        teta = fabs(a->phi(i,j,k))/(fabs(a->phi(i-1,j,k))+fabs(a->phi(i,j,k))) + 0.0001*p->DXN[IP]/(fabs(a->phi(i,j,k-1))+fabs(a->phi(i,j,k)));
        check=2;
        }
        
        if(p->flag4[IJp1K]==AIR)
        {
        teta = fabs(a->phi(i,j,k))/(fabs(a->phi(i,j+1,k))+fabs(a->phi(i,j,k))) + 0.0001*p->DYN[JP]/(fabs(a->phi(i,j+1,k))+fabs(a->phi(i,j,k)));
        check=1;
        }
        
        if(p->flag4[IJm1K]==AIR)
        {
        teta = fabs(a->phi(i,j,k))/(fabs(a->phi(i,j-1,k))+fabs(a->phi(i,j,k))) + 0.0001*p->DYN[JP]/(fabs(a->phi(i,j-1,k))+fabs(a->phi(i,j,k)));
        check=2;
        }
        
        if(p->flag4[IJKm1]==AIR)
        {
        teta = fabs(a->phi(i,j,k))/(fabs(a->phi(i,j,k-1))+fabs(a->phi(i,j,k))) + 0.0001*p->DZN[KP]/(fabs(a->phi(i,j,k-1))+fabs(a->phi(i,j,k)));
        check=2;
        }
        }
        
        if(p->flag4[IJKp1]==AIR)
        {
        teta = fabs(a->phi(i,j,k))/(fabs(a->phi(i,j,k+1))+fabs(a->phi(i,j,k))) + 0.0001*p->DZN[KP]/(fabs(a->phi(i,j,k+1))+fabs(a->phi(i,j,k)));
        check=1;
        }

    if(check==1)    
    wvel(i,j,k) -= alpha*p->dt*CPOR3*PORVAL3*(((1.0 - 1.0/teta)*a->press(i,j,k)-a->press(i,j,k))
    /(p->DZP[KP]*pd->roface(p,a,0,0,1)));
    
    if(check==2)    
    wvel(i,j,k) -= alpha*p->dt*CPOR3*PORVAL3*((a->press(i,j,k)-(1.0/teta)*a->press(i,j,k))
    /(p->DZP[KP]*pd->roface(p,a,0,0,1)));

        
    if(check==0)
    wvel(i,j,k) -= alpha*p->dt*CPOR3*PORVAL3*((a->press(i,j,k+1)-a->press(i,j,k))
	/(p->DZP[KP]*pd->roface(p,a,0,0,1)));
    }
}
 
void pjm_nse::rhs(lexer *p, fdm* a, ghostcell *pgc, field &u, field &v, field &w,double alpha)
{
	double uvel,vvel,wvel;
	
    count=0;
    FLUIDLOOP
    {
	a->rhsvec.V[count]=0.0;
    ++count;
    }
	
    pip=p->Y50;
    
    count=0;
    FLUIDLOOP
    {
    a->rhsvec.V[count] =  -(u(i,j,k)-u(i-1,j,k))/(alpha*p->dt*p->DXN[IP])
                          -(v(i,j,k)-v(i,j-1,k))/(alpha*p->dt*p->DYN[JP])
                          -(w(i,j,k)-w(i,j,k-1))/(alpha*p->dt*p->DZN[KP]);
                           
    ++count;
    }
    
    pip=0;
}
 
void pjm_nse::vel_setup(lexer *p, fdm* a, ghostcell *pgc, field &u, field &v, field &w,double alpha)
{
	pgc->start1(p,u,gcval_u);
	pgc->start2(p,v,gcval_v);
	pgc->start3(p,w,gcval_w);
}

void pjm_nse::upgrad(lexer*p,fdm* a, slice &eta, slice &eta_n)
{
}

void pjm_nse::vpgrad(lexer*p,fdm* a, slice &eta, slice &eta_n)
{
}

void pjm_nse::wpgrad(lexer*p,fdm* a, slice &eta, slice &eta_n)
{
}

void pjm_nse::ini(lexer*p,fdm* a, ghostcell *pgc){}


