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

#define HX (fabs(d->hx(i,j))>1.0e-20?d->hx(i,j):1.0e20)
#define HXP (fabs(0.5*(d->WL(i,j)+d->WL(i+1,j)))>1.0e-20?0.5*(d->WL(i,j)+d->WL(i+1,j)):1.0e20)
#define HY (fabs(d->hy(i,j))>1.0e-20?d->hy(i,j):1.0e20)

#include"nhflow_pjm.h"
#include"lexer.h"
#include"fdm_nhf.h" 
#include"ghostcell.h"
#include"solver.h"
#include"ioflow.h"
#include"nhflow_poisson.h"
#include"density_f.h"

 
nhflow_pjm::nhflow_pjm(lexer* p, fdm_nhf *d, ghostcell *pgc) : teta(1.0)
{
	pd = new density_f(p);
    
    ppois = new nhflow_poisson(p);
	
    gcval_press=540;  
}

nhflow_pjm::~nhflow_pjm()
{
}

void nhflow_pjm::start(fdm_nhf *d,lexer *p, poisson* ppois,solver* psolv, ghostcell* pgc, ioflow *pflow, double *U, double *V, double *W, double alpha)
{
    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<".";
			
    rhs(p,d,pgc,U,V,W,alpha);
    
    bedbc(p,d,pgc,U,V,W,alpha);
	
    ppois->start(p,d,d->P);
	
        starttime=pgc->timer();

    psolv->start(p,d,pgc,d->P,d->rhsvec,5);
	
        endtime=pgc->timer();
    
	pgc->start4(p,d->P,gcval_press);
	
	ucorr(p,d,uvel,alpha);
	vcorr(p,d,vvel,alpha);
	wcorr(p,d,wvel,alpha);

    p->poissoniter=p->solveriter;

	p->poissontime=endtime-starttime;

	if(p->mpirank==0 && p->count%p->P12==0)
	cout<<"piter: "<<p->solveriter<<"  ptime: "<<setprecision(3)<<p->poissontime<<endl;

    pgc->start4(p,d->test,1);
}

void nhflow_pjm::ucorr(lexer* p, fdm_nhf *d, double *U, double alpha)
{	
    if(p->D37==1)
	ULOOP
	uvel(i,j,k) -= alpha*p->dt*CPOR1*PORVAL1*(1.0/pd->roface(p,d,1,0,0))*((d->P(i+1,j,k)-d->P(i,j,k))/p->DXP[IP]
                + 0.25*(p->sigx[FIJK]+p->sigx[FIJKp1]+p->sigx[FIp1JK]+p->sigx[FIp1JKp1])*(0.5*(d->P(i,j,k+1)+d->P(i+1,j,k+1))-0.5*(d->P(i,j,k-1)+d->P(i+1,j,k-1)))/(p->DZP[KP]+p->DZP[KP1]));
    
    if(p->D37==2)
    ULOOP
    {       
     check=0;
    
        if(p->flag1[IJKp1]<0)
        check=1;        
    
    if(check==1)
    uvel(i,j,k) -= alpha*p->dt*CPOR1*PORVAL1*(1.0/pd->roface(p,d,1,0,0))*((d->P(i+1,j,k)-d->P(i,j,k))/p->DXP[IP]
                + 0.25*(p->sigx[FIJK]+p->sigx[FIJKp1]+p->sigx[FIp1JK]+p->sigx[FIp1JKp1])*(0.5*((1.0 - 1.0/teta)*(d->P(i,j,k)+d->P(i+1,j,k)))-0.5*(d->P(i,j,k-1)+d->P(i+1,j,k-1)))/(p->DZP[KP]+p->DZP[KP1]));
    
    if(check==0)
    uvel(i,j,k) -= alpha*p->dt*CPOR1*PORVAL1*(1.0/pd->roface(p,d,1,0,0))*((d->P(i+1,j,k)-d->P(i,j,k))/p->DXP[IP]
                + 0.25*(p->sigx[FIJK]+p->sigx[FIJKp1]+p->sigx[FIp1JK]+p->sigx[FIp1JKp1])*(0.5*(d->P(i,j,k+1)+d->P(i+1,j,k+1))-0.5*(d->P(i,j,k-1)+d->P(i+1,j,k-1)))/(p->DZP[KP]+p->DZP[KP1]));
    }
}

void nhflow_pjm::vcorr(lexer* p, fdm_nhf *d, double *V,double alpha)
{	
    if(p->D37==1)
    VLOOP
    vvel(i,j,k) -= alpha*p->dt*CPOR2*PORVAL2*(1.0/pd->roface(p,d,0,1,0))*((d->P(i,j+1,k)-d->P(i,j,k))/p->DYP[JP] 
                + 0.25*(p->sigy[FIJK]+p->sigy[FIJKp1]+p->sigy[FIJp1K]+p->sigy[FIJp1Kp1])*(0.5*(d->P(i,j,k+1)+d->P(i,j+1,k+1))-0.5*(d->P(i,j,k-1)+d->P(i,j+1,k-1)))/(p->DZP[KP]+p->DZP[KP1]));
                
                
    if(p->D37==2)
    VLOOP
    {       
     check=0;
    
        if(p->flag2[IJKp1]<0)
        check=1;        
    
    if(check==1)
    vvel(i,j,k) -= alpha*p->dt*CPOR2*PORVAL2*(1.0/pd->roface(p,d,0,1,0))*((d->P(i,j+1,k)-d->P(i,j,k))/p->DYP[JP] 
                + 0.25*(p->sigy[FIJK]+p->sigy[FIJKp1]+p->sigy[FIJp1K]+p->sigy[FIJp1Kp1])*(0.5*(1.0 - 1.0/teta)*(d->P(i,j,k)+d->P(i,j+1,k))-0.5*(d->P(i,j,k-1)+d->P(i,j+1,k-1)))/(p->DZP[KP]+p->DZP[KP1]));
    
    if(check==0)
    vvel(i,j,k) -= alpha*p->dt*CPOR2*PORVAL2*(1.0/pd->roface(p,d,0,1,0))*((d->P(i,j+1,k)-d->P(i,j,k))/p->DYP[JP] 
                + 0.25*(p->sigy[FIJK]+p->sigy[FIJKp1]+p->sigy[FIJp1K]+p->sigy[FIJp1Kp1])*(0.5*(d->P(i,j,k+1)+d->P(i,j+1,k+1))-0.5*(d->P(i,j,k-1)+d->P(i,j+1,k-1)))/(p->DZP[KP]+p->DZP[KP1]));
    }
}

void nhflow_pjm::wcorr(lexer* p, fdm_nhf *d, double *W, double alpha)
{
    if(p->D37==1)
    WLOOP 	
	wvel(i,j,k) -= alpha*p->dt*CPOR3*PORVAL3*((d->P(i,j,k+1)-d->P(i,j,k))/(p->DZP[KP]*pd->roface(p,d,0,0,1)))*p->sigz[IJ];
    
    
    if(p->D37==2)
	WLOOP
    {
    check=0;
    
        if(p->flag3[IJKp1]<0)
        check=1;

    if(check==1)    
    wvel(i,j,k) -= alpha*p->dt*CPOR3*PORVAL3*(((1.0 - 1.0/teta)*d->P(i,j,k)-d->P(i,j,k))/(p->DZP[KP]*pd->roface(p,d,0,0,1)))*p->sigz[IJ];
    
    if(check==0)
    wvel(i,j,k) -= alpha*p->dt*CPOR3*PORVAL3*((d->P(i,j,k+1)-d->P(i,j,k))/(p->DZP[KP]*pd->roface(p,d,0,0,1)))*p->sigz[IJ];
    }
}
 
void nhflow_pjm::rhs(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W, double alpha)
{
    NLOOP4
	d->rhsvec.V[n]=0.0;
	
    pip=p->Y50;

    n=0;
    LOOP
    {
    d->rhsvec.V[n] =      -  ((u(i,j,k)-u(i-1,j,k))/p->DXN[IP]
                            + 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*(0.5*(u(i,j,k+1)+u(i-1,j,k+1))-0.5*(u(i,j,k-1)+u(i-1,j,k-1)))/(p->DZP[KP]+p->DZP[KP1])
                            
                            + (v(i,j,k)-v(i,j-1,k))/p->DYN[JP] 
                            + 0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*(0.5*(v(i,j,k+1)+v(i,j-1,k+1))-0.5*(v(i,j,k-1)+v(i,j-1,k-1)))/(p->DZP[KP]+p->DZP[KP1])
                           
                            + p->sigz[IJ]*(w(i,j,k)-w(i,j,k-1))/p->DZN[KP])/(alpha*p->dt);
                           
    ++n;
    }
    pip=0;
}

void nhflow_pjm::bedbc(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W,double alpha)
{
    
}
 
void nhflow_pjm::vel_setup(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W,double alpha)
{

}

void nhflow_pjm::upgrad(lexer*p, fdm_nhf *d, slice &eta, slice &eta_n)
{
    if(p->D38==1 && p->A540==1)
    ULOOP
	d->F(i,j,k) -= PORVAL1*fabs(p->W22)*(p->A223*eta(i+1,j) + (1.0-p->A223)*eta_n(i+1,j) - p->A223*eta(i,j) - (1.0-p->A223)*eta_n(i,j))/p->DXP[IP];
    
    if(p->D38==1 && p->A540==2)
    ULOOP
	d->F(i,j,k) -= PORVAL1*fabs(p->W22)*(d->eta(i+1,j) - d->eta(i,j))/p->DXP[IP];
    
   
    if(p->D38==2 && p->A540==1)
    ULOOP
	d->F(i,j,k) -= PORVAL1*fabs(p->W22)*(1.0/HX)*
    
                    (0.5*(pow(eta(i+1,j),2.0) - pow(eta(i,j),2.0))/p->DXP[IP]
                    
                    + ((p->A223*eta(i+1,j) + (1.0-p->A223)*eta_n(i+1,j))*d->depth(i+1,j) - (p->A223*eta(i,j) + (1.0-p->A223)*eta_n(i,j))*d->depth(i,j))/p->DXP[IP]
                    
                    - 0.5*((p->A223*eta(i,j) + (1.0-p->A223)*eta_n(i,j)) + (p->A223*eta(i+1,j) + (1.0-p->A223)*eta_n(i+1,j)))*(d->depth(i+1,j)-d->depth(i,j))/p->DXP[IP]);
    
    if(p->D38==2 && p->A540==2)
    ULOOP
	d->F(i,j,k) -= PORVAL1*fabs(p->W22)*(1.0/HX)*
    
                    (0.5*(pow(d->eta(i+1,j),2.0) - pow(d->eta(i,j),2.0))/p->DXP[IP]
                    
                    + (d->eta(i+1,j)*d->depth(i+1,j) - d->eta(i,j)*d->depth(i,j))/p->DXP[IP]
                    
                    - 0.5*(d->eta(i,j) + d->eta(i+1,j))*(d->depth(i+1,j)-d->depth(i,j))/p->DXP[IP]);
    
    // fx = 1/2 g (eta^2 - 2* eta *z_b)
    // Sx = -g * eta * eta * Bx
}

void nhflow_pjm::vpgrad(lexer*p,fdm_nhf *d, slice &eta, slice &eta_n)
{
    if(p->D38==1 && p->A540==1)
    VLOOP
	d->G(i,j,k) -= PORVAL2*fabs(p->W22)*(p->A223*eta(i,j+1) + (1.0-p->A223)*eta_n(i,j+1) - p->A223*eta(i,j) - (1.0-p->A223)*eta_n(i,j))/p->DYP[JP];
    
    if(p->D38==1 && p->A540==2)
    VLOOP
	d->G(i,j,k) -= PORVAL2*fabs(p->W22)*(d->eta(i,j+1) - d->eta(i,j))/p->DYP[JP];
}

void nhflow_pjm::wpgrad(lexer*p, fdm_nhf *d, slice &eta, slice &eta_n)
{
}






