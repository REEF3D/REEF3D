/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"pjm_corr.h"
#include"lexer.h"
#include"fdm.h" 
#include"ghostcell.h"
#include"poisson.h"
#include"solver.h"
#include"momentum.h"
#include"ioflow.h"
#include"heat.h"
#include"concentration.h"
#include"density_f.h"
#include"density_comp.h"
#include"density_conc.h"
#include"density_heat.h"
#include"denisty_rheology.h"
#include"density_vof.h"
 
pjm_corr::pjm_corr(lexer* p, fdm *a, heat *&pheat, concentration *&pconc) : pcorr(p)
{
    if(p->F30>0 && p->H10==0 && p->W30==0 && p->W90==0)
	pd = new density_f(p);
	
	if(p->F30>0 && p->H10==0 && p->W30==1 && p->W90==0)
	pd = new density_comp(p);
	
	if(p->F30>0 && p->H10>0 && p->W90==0)
	pd = new density_heat(pheat);
	
	if(p->F30>0 && p->C10>0 && p->W90==0)
	pd = new density_conc(p,conc);
	
	if(p->F30>0 && p->H10==0 && p->W30==0 && p->W90>0)
	pd = new density_rheology(p);
    
    if(p->F80>0 && p->H10==0 && p->W30==0 && p->W90==0)
	pd = new density_vof(p);
    
    
    if(p->B76==0)
    gcval_press=40;  

    if(p->B76==1)
    gcval_press=41;

    if(p->B76==2)
    gcval_press=42;

    if(p->B76==3)
    gcval_press=43;
	
	if(p->B76==4) 
    gcval_press=44;
	
	if(p->B76==5) 
    gcval_press=45;
	
	gcval_u=7;
	gcval_v=8;
	gcval_w=9;
}

pjm_corr::~pjm_corr()
{
}

void pjm_corr::start(fdm* a,lexer*p, poisson* ppois,solver* psolv, ghostcell* pgc, momentum *pmom, ioflow *pflow, field& uvel, field& vvel, field& wvel, double alpha)
{
    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<".";
			
	vel_setup(p,a,pgc,uvel,vvel,wvel,alpha);	
    rhs(p,a,pgc,uvel,vvel,wvel,alpha);
    
    LOOP
    pcorr(i,j,k)=0.0;
    pgc->start4(p,pcorr,40);
	
    ppois->start(p,a,pcorr);
	
        starttime=pgc->timer();

    psolv->start(p,a,pgc,pcorr,a->xvec,a->rhsvec,5,gcval_press,p->N44);
	
        endtime=pgc->timer();
    
	pgc->start4(p,a->press,gcval_press);
    presscorr(p,a,pcorr,alpha);
    pgc->start4(p,pcorr,40);
	
	ucorr(a,p,uvel,alpha);
	vcorr(a,p,vvel,alpha);
	wcorr(a,p,wvel,alpha);

    p->poissoniter=p->solveriter;

	p->poissontime=endtime-starttime;

	if(p->mpirank==0 && innercounter==p->N50-1 && (p->count%p->P12==0))
	cout<<"piter: "<<p->solveriter<<"  ptime: "<<setprecision(3)<<p->poissontime<<endl;
    
}

void pjm_corr::ucorr(fdm* a, lexer* p, field& uvel,double alpha)
{	
	ULOOP
	uvel(i,j,k) -= alpha*p->dt*CPOR1*PORVAL1*((pcorr(i+1,j,k)-pcorr(i,j,k))
	/(p->DXP[IP]*pd->roface(p,a,1,0,0)));
}

void pjm_corr::vcorr(fdm* a, lexer* p, field& vvel,double alpha)
{	
	VLOOP
	vvel(i,j,k) -= alpha*p->dt*CPOR2*PORVAL2*((pcorr(i,j+1,k)-pcorr(i,j,k))
	/(p->DYP[JP]*pd->roface(p,a,0,1,0)));
}

void pjm_corr::wcorr(fdm* a, lexer* p, field& wvel,double alpha)
{	
	WLOOP
	wvel(i,j,k) -= alpha*p->dt*CPOR3*PORVAL3*((pcorr(i,j,k+1)-pcorr(i,j,k))
	/(p->DZP[KP]*pd->roface(p,a,0,0,1)));
}

void pjm_corr::presscorr(lexer* p, fdm* a, field& pcorr, double alpha)
{
    LOOP
    a->press(i,j,k) += pcorr(i,j,k);
}
 
void pjm_corr::rhs(lexer *p, fdm* a, ghostcell *pgc, field &u, field &v, field &w,double alpha)
{
    count=0;
	double uvel,vvel,wvel;
	
    NLOOP4
	a->rhsvec.V[n]=0.0;
	
    pip=p->Y50;

    LOOP
    {
    a->rhsvec.V[count] =  -(u(i,j,k)-u(i-1,j,k))/(alpha*p->dt*p->DXN[IP])
						   -(v(i,j,k)-v(i,j-1,k))/(alpha*p->dt*p->DYN[JP])
						   -(w(i,j,k)-w(i,j,k-1))/(alpha*p->dt*p->DZN[KP]);
                           
    ++count;
    }
    
    pip=0;
}
 

void pjm_corr::vel_setup(lexer *p, fdm* a, ghostcell *pgc, field &u, field &v, field &w,double alpha)
{
	pgc->start1(p,u,gcval_u);
	pgc->start2(p,v,gcval_v);
	pgc->start3(p,w,gcval_w);
	
	u.ggcpol(p);
	v.ggcpol(p);
	w.ggcpol(p);
}

void pjm_corr::pressure_norm(lexer*p, fdm* a, ghostcell* pgc)
{
    double sum=0.0;

    LOOP
    sum+=a->press(i,j,k);

    sum=pgc->globalsum(sum);

    sum/=double(p->cellnumtot);

    LOOP
    a->press(i,j,k)-=sum;
}

void pjm_corr::upgrad(lexer*p,fdm* a)
{
    ULOOP
    a->F(i,j,k)-=PORVAL1*(a->press(i+1,j,k)-a->press(i,j,k))/(p->DXP[IP]*pd->roface(p,a,1,0,0));
}

void pjm_corr::vpgrad(lexer*p,fdm* a)
{
    VLOOP
    a->G(i,j,k)-=PORVAL2*(a->press(i,j+1,k)-a->press(i,j,k))/(p->DYP[JP]*pd->roface(p,a,0,1,0));
}

void pjm_corr::wpgrad(lexer*p,fdm* a)
{
    WLOOP
    a->H(i,j,k)-=PORVAL3*(a->press(i,j,k+1)-a->press(i,j,k))/(p->DZP[KP]*pd->roface(p,a,0,0,1));
}

void pjm_corr::fillapu(lexer*p,fdm* a)
{
}

void pjm_corr::fillapv(lexer*p,fdm* a)
{
}

void pjm_corr::fillapw(lexer*p,fdm* a)
{
}

void pjm_corr::ptimesave(lexer *p, fdm *a, ghostcell *pgc)
{
}






