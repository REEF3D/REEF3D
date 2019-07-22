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

#include"pjm_sig.h"
#include"lexer.h"
#include"fdm.h" 
#include"ghostcell.h"
#include"poisson.h"
#include"solver.h"
#include"ioflow.h"
#include"heat.h"
#include"concentration.h"
#include"density_f.h"
#include"density_comp.h"
#include"density_conc.h"
#include"density_heat.h"
#include"density_rheology.h"
#include"density_vof.h"
 
pjm_sig::pjm_sig(lexer* p, fdm *a, heat *&pheat, concentration *&pconc)
{
    if(p->F30>0 && p->H10==0 && p->W30==0 && p->W90==0)
	pd = new density_f(p);
	
	if(p->F30>0 && p->H10==0 && p->W30==1 && p->W90==0)
	pd = new density_comp(p);
	
	if(p->F30>0 && p->H10>0 && p->W90==0)
	pd = new density_heat(p,pheat);
	
	if(p->F30>0 && p->C10>0 && p->W90==0)
	pd = new density_conc(p,pconc);
	
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

pjm_sig::~pjm_sig()
{
}

void pjm_sig::start(fdm* a,lexer*p, poisson* ppois,solver* psolv, ghostcell* pgc, momentum *pmom, ioflow *pflow, field& uvel, field& vvel, field& wvel, double alpha)
{
    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<".";
			
	vel_setup(p,a,pgc,uvel,vvel,wvel,alpha);	
    rhs(p,a,pgc,uvel,vvel,wvel,alpha);
	
    ppois->start(p,a,a->press);
	
        starttime=pgc->timer();

    psolv->start(p,a,pgc,a->press,a->xvec,a->rhsvec,5,gcval_press,p->N44);
	
        endtime=pgc->timer();
    
	pgc->start4(p,a->press,gcval_press);
	
	ucorr(p,a,uvel,alpha);
	vcorr(p,a,vvel,alpha);
	wcorr(p,a,wvel,alpha);
    
    p->poissoniter=p->solveriter;

	p->poissontime=endtime-starttime;

	if(p->mpirank==0 && innercounter==p->N50-1 && (p->count%p->P12==0))
	cout<<"piter: "<<p->solveriter<<"  ptime: "<<setprecision(3)<<p->poissontime<<endl;
}

void pjm_sig::ucorr(lexer* p, fdm* a, field& uvel,double alpha)
{	
	ULOOP
	uvel(i,j,k) -= alpha*p->dt*CPOR1*PORVAL1*((a->press(i+1,j,k)-a->press(i,j,k))
	/(p->DXP[IP]*pd->roface(p,a,1,0,0)) + p->sigx[FIJK]*(0.5*(a->press(i,j,k+1)+a->press(i+1,j,k+1))-0.5*(a->press(i,j,k-1)+a->press(i+1,j,k-1)))/(p->DZP[KP]+p->DZP[KP1]));
}

void pjm_sig::vcorr(lexer* p, fdm* a, field& vvel,double alpha)
{	 
    VLOOP
    vvel(i,j,k) -= alpha*p->dt*CPOR2*PORVAL2*(a->press(i,j+1,k)-a->press(i,j,k))
    /(p->DYP[JP]*(pd->roface(p,a,0,1,0)) + p->sigy[FIJK]*(0.5*(a->press(i,j,k+1)+a->press(i,j+1,k+1))-0.5*(a->press(i,j,k-1)+a->press(i,j+1,k-1)))/(p->DZP[KP]+p->DZP[KP1]));
}

void pjm_sig::wcorr(lexer* p, fdm* a, field& wvel,double alpha)
{	
	WLOOP
	wvel(i,j,k) -= alpha*p->dt*CPOR3*PORVAL3*((a->press(i,j,k+1)-a->press(i,j,k))
	/(p->DZP[KP]*pd->roface(p,a,0,0,1)))*p->sigz[IJ];
}
 
void pjm_sig::rhs(lexer *p, fdm* a, ghostcell *pgc, field &u, field &v, field &w,double alpha)
{
    NLOOP4
	a->rhsvec.V[n]=0.0;
	
    pip=p->Y50;

    count=0;
    LOOP
    {
    a->rhsvec.V[count] =  -(u(i,j,k)-u(i-1,j,k))/(alpha*p->dt*p->DXN[IP])  
                            + p->sigx[FIJK]*(0.5*(u(i,j,k+1)+u(i-1,j,k+1))-0.5*(u(i,j,k-1)+u(i-1,j,k-1)))/(p->DZP[KP]+p->DZP[KP1])
                            
						   -(v(i,j,k)-v(i,j-1,k))/(alpha*p->dt*p->DYN[JP])  
                           + p->sigy[FIJK]*(0.5*(v(i,j,k+1)+v(i,j-1,k+1))-0.5*(v(i,j,k-1)+v(i,j-1,k-1)))/(p->DZP[KP]+p->DZP[KP1])
                           
						   -(w(i,j,k)-w(i,j,k-1))/(alpha*p->dt*p->DZN[KP])*p->sigz[IJ];
    ++count;
    }
    pip=0;
}
 
void pjm_sig::vel_setup(lexer *p, fdm* a, ghostcell *pgc, field &u, field &v, field &w,double alpha)
{
	pgc->start1(p,u,gcval_u);
	pgc->start2(p,v,gcval_v);
	pgc->start3(p,w,gcval_w);
	
	u.ggcpol(p);
	v.ggcpol(p);
	w.ggcpol(p);
}

void pjm_sig::pressure_norm(lexer*p, fdm* a, ghostcell* pgc)
{
    double sum=0.0;

    LOOP
    sum+=a->press(i,j,k);

    sum=pgc->globalsum(sum);

    sum/=double(p->cellnumtot);

    LOOP
    a->press(i,j,k)-=sum;
}

void pjm_sig::upgrad(lexer*p,fdm* a)
{
}

void pjm_sig::vpgrad(lexer*p,fdm* a)
{
}

void pjm_sig::wpgrad(lexer*p,fdm* a)
{
}

void pjm_sig::fillapu(lexer*p,fdm* a)
{
}

void pjm_sig::fillapv(lexer*p,fdm* a)
{
}

void pjm_sig::fillapw(lexer*p,fdm* a)
{
}

void pjm_sig::ptimesave(lexer *p, fdm *a, ghostcell *pgc)
{
}






