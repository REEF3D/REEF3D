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

#define HX (fabs(a->hx(i,j))>1.0e-20?a->hx(i,j):1.0e20)
#define HXP (fabs(0.5*(a->WL(i,j)+a->WL(i+1,j)))>1.0e-20?0.5*(a->WL(i,j)+a->WL(i+1,j)):1.0e20)
#define HY (fabs(a->hy(i,j))>1.0e-20?a->hy(i,j):1.0e20)

#include"nhflow_pjmcorr.h"
#include"lexer.h"
#include"fdm_nhf.h" 
#include"ghostcell.h"
#include"nhflow_poisson.h"
#include"solver.h"
#include"ioflow.h"
#include"heat.h"
#include"concentration.h"
#include"density_f.h"
#include"density_comp.h"
#include"density_conc.h"
#include"density_heat.h"
#include"density_vof.h"
 
nhflow_pjmcorr::nhflow_pjmcorr(lexer* p, fdm_nhf *d, ghostcell *pgc) : teta(0.5)
{
	pd = new density_f(p);
	
    gcval_press=540;  
}

nhflow_pjmcorr::~nhflow_pjmcorr()
{
}

void nhflow_pjmcorr::start(lexer*p, fdm_nhf *d, solver* psolv, ghostcell* pgc, ioflow *pflow, 
                            double *U, double *V, double *W, double alpha)
{
    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<".";
			
    rhs(p,a,pgc,uvel,vvel,wvel,alpha);
    
    bedbc(p,a,pgc,uvel,vvel,wvel,alpha);
	
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

	if(p->mpirank==0 && p->count%p->P12==0)
	cout<<"piter: "<<p->solveriter<<"  ptime: "<<setprecision(3)<<p->poissontime<<endl;

    pgc->start4(p,a->test,1);
}

void nhflow_pjmcorr::ucorr(lexer* p, fdm_nhf *d, double *U, double alpha)
{	
    if(p->D37==1)
	LOOP
	U[IJK] -= alpha*p->dt*CPOR1*PORVAL1*(1.0/pd->roface(p,a,1,0,0))*((a->press(i+1,j,k)-a->press(i,j,k))/p->DXP[IP]
                + 0.25*(p->sigx[FIJK]+p->sigx[FIJKp1]+p->sigx[FIp1JK]+p->sigx[FIp1JKp1])*(0.5*(a->press(i,j,k+1)+a->press(i+1,j,k+1))-0.5*(a->press(i,j,k-1)+a->press(i+1,j,k-1)))/(p->DZP[KP]+p->DZP[KP1]));
    
    if(p->D37==2)
    LOOP
    {       
     check=0;
    
        if(p->flag1[IJKp1]<0)
        check=1;        
    
    if(check==1)
    U[IJK] -= alpha*p->dt*CPOR1*PORVAL1*(1.0/pd->roface(p,a,1,0,0))*((a->press(i+1,j,k)-a->press(i,j,k))/p->DXP[IP]
                + 0.25*(p->sigx[FIJK]+p->sigx[FIJKp1]+p->sigx[FIp1JK]+p->sigx[FIp1JKp1])*(0.5*((1.0 - 1.0/teta)*(a->press(i,j,k)+a->press(i+1,j,k)))-0.5*(a->press(i,j,k-1)+a->press(i+1,j,k-1)))/(p->DZP[KP]+p->DZP[KP1]));
    
    if(check==0)
    U[IJK] -= alpha*p->dt*CPOR1*PORVAL1*(1.0/pd->roface(p,a,1,0,0))*((a->press(i+1,j,k)-a->press(i,j,k))/p->DXP[IP]
                + 0.25*(p->sigx[FIJK]+p->sigx[FIJKp1]+p->sigx[FIp1JK]+p->sigx[FIp1JKp1])*(0.5*(a->press(i,j,k+1)+a->press(i+1,j,k+1))-0.5*(a->press(i,j,k-1)+a->press(i+1,j,k-1)))/(p->DZP[KP]+p->DZP[KP1]));
    }
}

void nhflow_pjmcorr::vcorr(lexer* p, fdm_nhf *d, double *V, double alpha)
{	
    if(p->D37==1)
    LOOP
    V[IJK] -= alpha*p->dt*CPOR2*PORVAL2*(1.0/pd->roface(p,a,0,1,0))*((a->press(i,j+1,k)-a->press(i,j,k))/p->DYP[JP] 
                + 0.25*(p->sigy[FIJK]+p->sigy[FIJKp1]+p->sigy[FIJp1K]+p->sigy[FIJp1Kp1])*(0.5*(a->press(i,j,k+1)+a->press(i,j+1,k+1))-0.5*(a->press(i,j,k-1)+a->press(i,j+1,k-1)))/(p->DZP[KP]+p->DZP[KP1]));
                
                
    if(p->D37==2)
    LOOP
    {       
     check=0;
    
        if(p->flag2[IJKp1]<0)
        check=1;        
    
    if(check==1)
    V[IJK] -= alpha*p->dt*CPOR2*PORVAL2*(1.0/pd->roface(p,a,0,1,0))*((a->press(i,j+1,k)-a->press(i,j,k))/p->DYP[JP] 
                + 0.25*(p->sigy[FIJK]+p->sigy[FIJKp1]+p->sigy[FIJp1K]+p->sigy[FIJp1Kp1])*(0.5*(1.0 - 1.0/teta)*(a->press(i,j,k)+a->press(i,j+1,k))-0.5*(a->press(i,j,k-1)+a->press(i,j+1,k-1)))/(p->DZP[KP]+p->DZP[KP1]));
    
    if(check==0)
    V[IJK] -= alpha*p->dt*CPOR2*PORVAL2*(1.0/pd->roface(p,a,0,1,0))*((a->press(i,j+1,k)-a->press(i,j,k))/p->DYP[JP] 
                + 0.25*(p->sigy[FIJK]+p->sigy[FIJKp1]+p->sigy[FIJp1K]+p->sigy[FIJp1Kp1])*(0.5*(a->press(i,j,k+1)+a->press(i,j+1,k+1))-0.5*(a->press(i,j,k-1)+a->press(i,j+1,k-1)))/(p->DZP[KP]+p->DZP[KP1]));
    }
}

void nhflow_pjmcorr::wcorr(lexer* p, fdm_nhf *d, double *W, double alpha)
{
    if(p->D37==1)
    LOOP 	
	W[IJK] -= alpha*p->dt*CPOR3*PORVAL3*((a->press(i,j,k+1)-a->press(i,j,k))/(p->DZP[KP]*pd->roface(p,a,0,0,1)))*p->sigz[IJ];
    
    
    if(p->D37==2)
	LOOP
    {
    check=0;
    
        if(p->flag3[IJKp1]<0)
        check=1;

    if(check==1)    
    W[IJK] -= alpha*p->dt*CPOR3*PORVAL3*(((1.0 - 1.0/teta)*a->press(i,j,k)-a->press(i,j,k))/(p->DZP[KP]*pd->roface(p,a,0,0,1)))*p->sigz[IJ];
    
    if(check==0)
    W[IJK] -= alpha*p->dt*CPOR3*PORVAL3*((a->press(i,j,k+1)-a->press(i,j,k))/(p->DZP[KP]*pd->roface(p,a,0,0,1)))*p->sigz[IJ];
    }
}
 
void nhflow_pjmcorr::rhs(lexer *p, fdm_nhf *d, ghostcell *pgc, field &u, field &v, field &w,double alpha)
{
    NLOOP4
	a->rhsvec.V[n]=0.0;
	
    pip=p->Y50;

    n=0;
    LOOP
    {
    a->rhsvec.V[n] =      -  ((u(i,j,k)-u(i-1,j,k))/p->DXN[IP]
                            + 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*(0.5*(u(i,j,k+1)+u(i-1,j,k+1))-0.5*(u(i,j,k-1)+u(i-1,j,k-1)))/(p->DZP[KP]+p->DZP[KP1])
                            
                            + (v(i,j,k)-v(i,j-1,k))/p->DYN[JP] 
                            + 0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*(0.5*(v(i,j,k+1)+v(i,j-1,k+1))-0.5*(v(i,j,k-1)+v(i,j-1,k-1)))/(p->DZP[KP]+p->DZP[KP1])
                           
                            + p->sigz[IJ]*(w(i,j,k)-w(i,j,k-1))/p->DZN[KP])/(alpha*p->dt);
                           
    ++n;
    }
    pip=0;
}

void nhflow_pjmcorr::bedbc(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W, double alpha)
{
    
}
 
void nhflow_pjmcorr::vel_setup(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W, double alpha)
{

}

void nhflow_pjmcorr::upgrad(lexer*p,fdm_nhf *d, slice &eta, slice &eta_n)
{
    if(p->D38==1 && p->A540==1)
    LOOP
	d->F[IJK] -= PORVAL1*fabs(p->W22)*(p->A223*eta(i+1,j) + (1.0-p->A223)*eta_n(i+1,j) - p->A223*eta(i,j) - (1.0-p->A223)*eta_n(i,j))/p->DXP[IP];
    
    if(p->D38==1 && p->A540==2)
    LOOP
	d->F[IJK] -= PORVAL1*fabs(p->W22)*(a->eta(i+1,j) - a->eta(i,j))/p->DXP[IP];
    
   
    if(p->D38==2 && p->A540==1)
    LOOP
	d->F[IJK] -= PORVAL1*fabs(p->W22)*(1.0/HX)*
    
                    (0.5*(pow(eta(i+1,j),2.0) - pow(eta(i,j),2.0))/p->DXP[IP]
                    
                    + ((p->A223*eta(i+1,j) + (1.0-p->A223)*eta_n(i+1,j))*a->depth(i+1,j) - (p->A223*eta(i,j) + (1.0-p->A223)*eta_n(i,j))*a->depth(i,j))/p->DXP[IP]
                    
                    - 0.5*((p->A223*eta(i,j) + (1.0-p->A223)*eta_n(i,j)) + (p->A223*eta(i+1,j) + (1.0-p->A223)*eta_n(i+1,j)))*(a->depth(i+1,j)-a->depth(i,j))/p->DXP[IP]);
    
    if(p->D38==2 && p->A540==2)
    LOOP
	d->F[IJK] -= PORVAL1*fabs(p->W22)*(1.0/HX)*
    
                    (0.5*(pow(a->eta(i+1,j),2.0) - pow(a->eta(i,j),2.0))/p->DXP[IP]
                    
                    + (a->eta(i+1,j)*a->depth(i+1,j) - a->eta(i,j)*a->depth(i,j))/p->DXP[IP]
                    
                    - 0.5*(a->eta(i,j) + a->eta(i+1,j))*(a->depth(i+1,j)-a->depth(i,j))/p->DXP[IP]);
    
    // fx = 1/2 g (eta^2 - 2* eta *z_b)
    // Sx = -g * eta * eta * Bx
}

void nhflow_pjmcorr::vpgrad(lexer*p,fdm_nhf *d, slice &eta, slice &eta_n)
{
    if(p->D38==1 && p->A540==1)
    LOOP
	d->G[IJK] -= PORVAL2*fabs(p->W22)*(p->A223*eta(i,j+1) + (1.0-p->A223)*eta_n(i,j+1) - p->A223*eta(i,j) - (1.0-p->A223)*eta_n(i,j))/p->DYP[JP];
    
    if(p->D38==1 && p->A540==2)
    LOOP
	d->G[IJK] -= PORVAL2*fabs(p->W22)*(a->eta(i,j+1) - a->eta(i,j))/p->DYP[JP];
}

void nhflow_pjmcorr::wpgrad(lexer*p,fdm_nhf *d, slice &eta, slice &eta_n)
{
}






