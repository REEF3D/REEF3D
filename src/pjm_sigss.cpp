/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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

#include"pjm_sigss.h"
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
#include"density_vof.h"
#include"hypre_struct.h"
#include"hypre_sstruct_fnpf.h"
 
pjm_sigss::pjm_sigss(lexer* p, fdm *a, ghostcell *pgc, heat *&pheat, concentration *&pconc)
{
    if((p->F80==0||p->A10==5||p->A10==55) && p->H10==0 && p->W30==0)
	pd = new density_f(p);
	
	if(p->F80==0 && p->H10==0 && p->W30==1)
	pd = new density_comp(p);
	
	if(p->F80==0 && p->H10>0)
	pd = new density_heat(p,pheat);
	
	if(p->F80==0 && p->C10>0)
	pd = new density_conc(p,pconc);

    gcval_press=540;  

	gcval_u=7;
	gcval_v=8;
	gcval_w=9;

    vecsize=p->knox*p->knoy*(p->knoz+1); 
    
    p->Darray(M,vecsize*15);
    p->Darray(x,vecsize);
    p->Darray(rhs,vecsize);
}

pjm_sigss::~pjm_sigss()
{
}

void pjm_sigss::start(fdm* a,lexer*p, poisson* ppois,solver* psolv, ghostcell* pgc, ioflow *pflow, field& uvel, field& vvel, field& wvel, double alpha)
{
    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<".";
			
	vel_setup(p,a,pgc,uvel,vvel,wvel,alpha);	
    rhscalc(p,a,pgc,uvel,vvel,wvel,alpha);
    
    if(p->j_dir==0)
    poisson2D(p,a,a->press);
    
    if(p->j_dir==1)
    poisson3D(p,a,a->press);
    
    fillvec(p,a,pgc);
	
        starttime=pgc->timer();

    psolv->startM(p,a,pgc,x,rhs,M,5);

        endtime=pgc->timer();
        
    fillvec_back(p,a,pgc);
    
	pgc->start4(p,a->press,gcval_press);
    
	ucorr(p,a,uvel,alpha);
	vcorr(p,a,vvel,alpha);
	wcorr(p,a,wvel,alpha);

    p->poissoniter=p->solveriter;

	p->poissontime=endtime-starttime;

	if(p->mpirank==0 && p->count%p->P12==0)
	cout<<"piter: "<<p->solveriter<<"  ptime: "<<setprecision(3)<<p->poissontime<<endl;
}

void pjm_sigss::ucorr(lexer* p, fdm* a, field& uvel,double alpha)
{	
	ULOOP
	uvel(i,j,k) -= alpha*p->dt*CPOR1*PORVAL1*(1.0/pd->roface(p,a,1,0,0))*((a->press(i+1,j,k)-a->press(i,j,k))/p->DXP[IP]
                + 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*(0.5*(a->press(i,j,k+1)+a->press(i+1,j,k+1))-0.5*(a->press(i,j,k-1)+a->press(i+1,j,k-1)))/(p->DZP[KP]+p->DZP[KP1]));
}

void pjm_sigss::vcorr(lexer* p, fdm* a, field& vvel,double alpha)
{	 
    VLOOP
    vvel(i,j,k) -= alpha*p->dt*CPOR2*PORVAL2*(1.0/pd->roface(p,a,0,1,0))*((a->press(i,j+1,k)-a->press(i,j,k))/p->DYP[JP] 
                + 0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*(0.5*(a->press(i,j,k+1)+a->press(i,j+1,k+1))-0.5*(a->press(i,j,k-1)+a->press(i,j+1,k-1)))/(p->DZP[KP]+p->DZP[KP1]));
}

void pjm_sigss::wcorr(lexer* p, fdm* a, field& wvel,double alpha)
{
    WLOOP 	
	wvel(i,j,k) -= alpha*p->dt*CPOR3*PORVAL3*((a->press(i,j,k+1)-a->press(i,j,k))/(p->DZP[KP]*pd->roface(p,a,0,0,1)))*p->sigz[IJ];
}
 
void pjm_sigss::rhscalc(lexer *p, fdm* a, ghostcell *pgc, field &u, field &v, field &w,double alpha)
{
    NLOOP4
	rhs[n]=0.0;
	
    pip=p->Y50;

    n=0;
    KJILOOP
    {
    PCHECK
    rhs[n] =        - ((u(i,j,k)-u(i-1,j,k))/p->DXN[IP]
                     + 0.25*(p->sigx[FIJK]+p->sigx[FIJKp1]+p->sigx[FIp1JK]+p->sigx[FIp1JKp1])*(0.5*(u(i,j,k+1)+u(i-1,j,k+1))-0.5*(u(i,j,k-1)+u(i-1,j,k-1)))/(p->DZP[KP]+p->DZP[KP1])
                            
                     + (v(i,j,k)-v(i,j-1,k))/p->DYN[JP] 
                     + 0.25*(p->sigy[FIJK]+p->sigy[FIJKp1]+p->sigy[FIJp1K]+p->sigy[FIJp1Kp1])*(0.5*(v(i,j,k+1)+v(i,j-1,k+1))-0.5*(v(i,j,k-1)+v(i,j-1,k-1)))/(p->DZP[KP]+p->DZP[KP1])
                           
                     + p->sigz[IJ]*(w(i,j,k)-w(i,j,k-1))/p->DZN[KP] )/(alpha*p->dt);
                           
    SCHECK
    rhs[n] = 0.0;
                                        
    ++n;
    }
    pip=0;
}
 
void pjm_sigss::vel_setup(lexer *p, fdm* a, ghostcell *pgc, field &u, field &v, field &w,double alpha)
{
	pgc->start1(p,u,gcval_u);
	pgc->start2(p,v,gcval_v);
	pgc->start3(p,w,gcval_w);
}

void pjm_sigss::upgrad(lexer*p,fdm* a, slice &eta, slice &eta_n)
{
    ULOOP
	a->F(i,j,k) -= PORVAL1*fabs(p->W22)*(p->A223*eta(i+1,j) + (1.0-p->A223)*eta_n(i+1,j) - p->A223*eta(i,j) - (1.0-p->A223)*eta_n(i,j))/p->DXP[IP];
}

void pjm_sigss::vpgrad(lexer*p,fdm* a, slice &eta, slice &eta_n)
{
    VLOOP
	a->G(i,j,k) -= PORVAL2*fabs(p->W22)*(p->A223*eta(i,j+1) + (1.0-p->A223)*eta_n(i,j+1) - p->A223*eta(i,j) - (1.0-p->A223)*eta_n(i,j))/p->DYP[JP];
}

void pjm_sigss::wpgrad(lexer*p,fdm* a, slice &eta, slice &eta_n)
{
}

void pjm_sigss::fillvec(lexer *p, fdm *a, ghostcell *pgc)
{
    n=0;
    KJILOOP
    {
    PCHECK
    x[n] = a->press(i,j,k);
    
    SCHECK
    x[n] = 0.0;
    
    ++n;
    }
}

void pjm_sigss::fillvec_back(lexer *p, fdm *a, ghostcell *pgc)
{
    n=0;
    KJILOOP
    {
    PCHECK
    a->press(i,j,k) = x[n];
    
    ++n;
    }
}






