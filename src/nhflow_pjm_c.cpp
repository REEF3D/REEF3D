/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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
#define WLVL (fabs(d->WL(i,j))>1.0e-20?d->WL(i,j):1.0e20)

#include"nhflow_pjm_c.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"nhflow_poisson_c.h"
#include"solver.h"
#include"ioflow.h"
#include"density_f.h"
#include"patchBC_interface.h"

nhflow_pjm_c::nhflow_pjm_c(lexer* p, fdm_nhf *d, ghostcell *pgc, patchBC_interface *ppBC) : teta(0.5)
{
    pBC = ppBC;
    
	pd = new density_f(p);

    ppois = new nhflow_poisson_c(p);

    gcval_press=540;
    
    if(p->D33==0)
    solver_id = 8;
    
    if(p->D33==1)
    solver_id = 9;
}

nhflow_pjm_c::~nhflow_pjm_c()
{
}

void nhflow_pjm_c::start(lexer *p, fdm_nhf *d, solver* psolv, ghostcell* pgc, ioflow *pflow,
                        double *U, double *V, double *W, double alpha)
{
    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<".";

    rhs(p,d,pgc,U,V,W,alpha);

    ppois->start(p,d,d->P);

        starttime=pgc->timer();

    psolv->startF(p,pgc,d->P,d->rhsvec,d->M,7);
    //psolv->start(p,a,pgc,a->press,a->rhsvec,5);

        endtime=pgc->timer();

	pgc->start4V(p,d->P,gcval_press);
    
	ucorr(p,d,U,alpha);
	//vcorr(p,d,V,alpha);
	wcorr(p,d,W,alpha);

    p->poissoniter=p->solveriter;

	p->poissontime=endtime-starttime;

	if(p->mpirank==0 && p->count%p->P12==0)
	cout<<"piter: "<<p->solveriter<<"  ptime: "<<setprecision(3)<<p->poissontime<<endl;
}

void nhflow_pjm_c::ucorr(lexer* p, fdm_nhf *d, double *U, double alpha)
{
    if(p->D37==1)
	ULOOP
	U[IJK] -= alpha*p->dt*CPORNH*PORVALNH*(1.0/p->W1)*((d->P[Ip1JK]-d->P[IJK])/p->DXP[IP]
    
                + 0.25*(p->sigx[FIJK]+p->sigx[FIJKp1]+p->sigx[FIp1JK]+p->sigx[FIp1JKp1])*
                
                (0.5*(d->P[IJKp1]+d->P[Ip1JKp1])-0.5*(d->P[IJKm1]+d->P[Ip1JKm1]))/(p->DZP[KP]+p->DZP[KP1]));
    
    if(p->D37>=2)
    ULOOP
    {       
     check=0;
    
        if(p->flag1[IJKp1]<0)
        check=1;        
    
    if(check==1)
    U[IJK] -= alpha*p->dt*CPORNH*PORVALNH*(1.0/p->W1)*((d->P[Ip1JK]-d->P[IJK])/p->DXP[IP]
                + 0.25*(p->sigx[FIJK]+p->sigx[FIJKp1]+p->sigx[FIp1JK]+p->sigx[FIp1JKp1])*(0.5*((1.0 - 1.0/teta)*(d->P[IJK]+d->P[Ip1JK]))-0.5*(d->P[IJKm1]+d->P[Ip1JKm1]))/(p->DZP[KP]+p->DZP[KP1]));
    
    if(check==0)
    U[IJK] -= alpha*p->dt*CPORNH*PORVALNH*(1.0/p->W1)*((d->P[Ip1JK]-d->P[IJK])/p->DXP[IP]
                + 0.25*(p->sigx[FIJK]+p->sigx[FIJKp1]+p->sigx[FIp1JK]+p->sigx[FIp1JKp1])*(0.5*(d->P[IJKp1]+d->P[Ip1JKp1])-0.5*(d->P[IJKm1]+d->P[Ip1JKm1]))/(p->DZP[KP]+p->DZP[KP1]));
    }
}

void nhflow_pjm_c::vcorr(lexer* p, fdm_nhf *d, double *V, double alpha)
{
    //if(p->D37==1)
    VLOOP
    V[IJK] -= alpha*p->dt*CPORNH*PORVALNH*(1.0/p->W1)*((d->P[IJp1K]-d->P[IJK])/p->DYP[JP] 
                + 0.25*(p->sigy[FIJK]+p->sigy[FIJKp1]+p->sigy[FIJp1K]+p->sigy[FIJp1Kp1])*(0.5*(d->P[IJKp1]+d->P[IJp1Kp1])-0.5*(d->P[IJKm1]+d->P[IJp1Km1]))/(p->DZP[KP]+p->DZP[KP1]));
            
                
   /* if(p->D37>=2)
    VLOOP
    {       
     check=0;
    
        if(p->flag2[IJKp1]<0)
        check=1;        
    
    if(check==1)
    vvel(i,j,k) -= alpha*p->dt*CPOR2*PORVAL2*(1.0/pd->roface(p,a,0,1,0))*((d->P[IJp1K]-d->P[IJK])/p->DYP[JP] 
                + 0.25*(p->sigy[FIJK]+p->sigy[FIJKp1]+p->sigy[FIJp1K]+p->sigy[FIJp1Kp1])*(0.5*(1.0 - 1.0/teta)*(d->P[IJK]+d->P[IJp1K])-0.5*(a->press(i,j,k-1)+a->press(i,j+1,k-1)))/(p->DZP[KP]+p->DZP[KP1]));
    
    if(check==0)
    vvel(i,j,k) -= alpha*p->dt*CPOR2*PORVAL2*(1.0/pd->roface(p,a,0,1,0))*((d->P[IJp1K]-d->P[IJK])/p->DYP[JP] 
                + 0.25*(p->sigy[FIJK]+p->sigy[FIJKp1]+p->sigy[FIJp1K]+p->sigy[FIJp1Kp1])*(0.5*(d->P[IJKp1]+a->press(i,j+1,k+1))-0.5*(a->press(i,j,k-1)+a->press(i,j+1,k-1)))/(p->DZP[KP]+p->DZP[KP1]));
    }*/
}

void nhflow_pjm_c::wcorr(lexer* p, fdm_nhf *d, double *W, double alpha)
{
    
    if(p->D37==1)
    WLOOP
	W[IJK] -= alpha*p->dt*CPORNH*PORVALNH*((d->P[IJKp1]-d->P[IJK])/(p->DZP[KP]*p->W1))*p->sigz[IJ];
    
    
    if(p->D37>=2)
	WLOOP
    {
    check=0;
    
        if(p->flag3[IJKp1]<0)
        check=1;

    if(check==1)    
    W[IJK] -= alpha*p->dt*CPORNH*PORVALNH*(((1.0 - 1.0/teta)*d->P[IJK] - d->P[IJK])/(p->DZP[KP]*p->W1))*p->sigz[IJ];
    
    if(check==0)
    W[IJK] -= alpha*p->dt*CPORNH*PORVALNH*((d->P[IJKp1]-d->P[IJK])/(p->DZP[KP]*p->W1))*p->sigz[IJ];
    }
    
    /*
    if(p->D37>=2)
	WLOOP
    {
    check=0;
    
        if(p->flag3[IJKp1]<0)
        check=1;

    if(check==1)    
    W[IJK] -= alpha*p->dt*CPOR3*PORVAL3*(((1.0 - 1.0/teta)*d->P[IJK]-d->P[IJK])/(p->DZP[KP]*pd->roface(p,a,0,0,1)))*p->sigz[IJ];
    
    if(check==0)
    W[IJK] -= alpha*p->dt*CPOR3*PORVAL3*((d->P[IJKp1]-d->P[IJK])/(p->DZP[KP]*pd->roface(p,a,0,0,1)))*p->sigz[IJ];
    }*/
}

void nhflow_pjm_c::rhs(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W, double alpha)
{
    n=0;
    FLOOP
    {
	d->rhsvec.V[n]=0.0;
    ++n;
    }

    n=0;
    LOOP
    {
    d->rhsvec.V[n] =      -  ((U[IJK]-U[Im1JK])/p->DXN[IP]
                            + 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*(0.5*(U[IJKp1]+U[Im1JKp1])-0.5*(U[IJKm1]+U[Im1JKm1]))/(p->DZP[KP]+p->DZP[KP1])
                            
                            + (V[IJK]-V[IJm1K])/p->DYN[JP] 
                            + 0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*(0.5*(V[IJKp1]+V[IJm1Kp1])-0.5*(V[IJKm1]+V[IJm1Km1]))/(p->DZP[KP]+p->DZP[KP1])
                           
                            + p->sigz[IJ]*(W[IJK]-W[IJKm1])/p->DZN[KP])/(alpha*p->dt);
                           
    ++n;
    }
    
    n=0;
    LOOP
    {
    d->test[IJK] =      -  ((U[IJK]-U[Im1JK])/p->DXN[IP]
                            + 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*(0.5*(U[IJKp1]+U[Im1JKp1])-0.5*(U[IJKm1]+U[Im1JKm1]))/(p->DZP[KP]+p->DZP[KP1])
                            
                            + (V[IJK]-V[IJm1K])/p->DYN[JP] 
                            + 0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*(0.5*(V[IJKp1]+V[IJm1Kp1])-0.5*(V[IJKm1]+V[IJm1Km1]))/(p->DZP[KP]+p->DZP[KP1])
                           
                            + p->sigz[IJ]*(W[IJK]-W[IJKm1])/p->DZN[KP])/(alpha*p->dt);
                           
    ++n;
    }
  
}

void nhflow_pjm_c::bedbc(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W,double alpha)
{
}

void nhflow_pjm_c::vel_setup(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W,double alpha)
{
}

void nhflow_pjm_c::upgrad(lexer*p, fdm_nhf *d, slice &eta, slice &eta_n)
{
    if(p->A521==1 && p->A540==1)
    ULOOP
    WETDRY
    d->F[IJK] -= PORVALNH*fabs(p->W22)*
                (p->A523*eta(i+1,j) + (1.0-p->A523)*eta_n(i+1,j) - p->A523*eta(i,j) - (1.0-p->A523)*eta_n(i,j))/(p->DXP[IP]);

    if(p->A521==1 && p->A540==2)
    ULOOP
    WETDRY
	d->F[IJK] -= PORVALNH*fabs(p->W22)*(d->eta(i+1,j) - d->eta(i,j))/p->DXP[IP];
               
    if(p->A521==2 && p->A540==1)
    ULOOP
    WETDRY
    if(p->wet[Ip1J]==1 || p->wet[Im1J]==1)
    {
        if(p->wet[Ip1J]==1 && p->wet[Im1J]==1)
        {
        dfdx_plus = (eta(i+1,j)-eta(i,j))/p->DXP[IP];
        dfdx_min  = (eta(i,j)-eta(i-1,j))/p->DXP[IM1];
        
        detadx = limiter(dfdx_plus,dfdx_min);
        
        dfdx_plus = (eta_n(i+1,j)-eta_n(i,j))/p->DXP[IP];
        dfdx_min  = (eta_n(i,j)-eta_n(i-1,j))/p->DXP[IM1];
        
        detadx_n = limiter(dfdx_plus,dfdx_min);
            
        d->F[IJK] -= PORVALNH*fabs(p->W22)*
                    (p->A523*detadx + (1.0-p->A523)*detadx_n);
        }
        
        if(p->wet[Ip1J]==0 && p->wet[Im1J]==1)
        {
        detadx = (eta(i,j)-eta(i-1,j))/p->DXP[IM1];
        detadx_n = (eta_n(i,j)-eta_n(i-1,j))/p->DXP[IM1];
            
        d->F[IJK] -= PORVALNH*fabs(p->W22)*
                    (p->A523*detadx + (1.0-p->A523)*detadx_n);
        }
        
        if(p->wet[Ip1J]==0 && p->wet[Im1J]==1)
        {
        detadx = (eta(i+1,j)-eta(i,j))/p->DXP[IP];
        detadx_n = (eta_n(i+1,j)-eta_n(i,j))/p->DXP[IP];
            
        d->F[IJK] -= PORVALNH*fabs(p->W22)*
                    (p->A523*detadx + (1.0-p->A523)*detadx_n);
        }
    }
}

void nhflow_pjm_c::vpgrad(lexer*p,fdm_nhf *d, slice &eta, slice &eta_n)
{
    if(p->A521==1 && p->A540==1)
    VLOOP
    WETDRY
	d->G[IJK] -= PORVALNH*fabs(p->W22)*
                 (p->A523*eta(i,j+1) + (1.0-p->A523)*eta_n(i,j+1) - p->A523*eta(i,j) - (1.0-p->A523)*eta_n(i,j))/(p->DYP[JP]);
    
    if(p->A521==1 && p->A540==2)
    VLOOP
    WETDRY
	d->G[IJK] -= PORVALNH*fabs(p->W22)*(d->eta(i,j+1) - d->eta(i,j))/p->DYP[JP];
    
    if(p->A521==2 && p->A540==1)
    VLOOP
    WETDRY
    if(p->wet[IJp1]==1 || p->wet[IJm1]==1)
    {
        if(p->wet[IJm1]==1 && p->wet[IJm1]==1)
        {
        dfdy_plus = (eta(i,j+1)-eta(i,j))/p->DYP[JP];
        dfdy_min  = (eta(i,j)-eta(i,j-1))/p->DYP[JM1];
        
        detady = limiter(dfdy_plus,dfdy_min);
        
        dfdy_plus = (eta_n(i,j+1)-eta_n(i,j))/p->DYP[JP];
        dfdy_min  = (eta_n(i,j)-eta_n(i,j-1))/p->DYP[JM1];
        
        detady_n = limiter(dfdy_plus,dfdy_min);
            
        d->G[IJK] -= PORVALNH*fabs(p->W22)*
                    (p->A523*detady + (1.0-p->A523)*detady_n);
        }
        
        if(p->wet[IJp1]==0 && p->wet[IJm1]==1)
        {
        detady = (eta(i,j)-eta(i,j-1))/p->DYP[JM1];
        detady_n = (eta_n(i,j)-eta_n(i,j-1))/p->DYP[JM1];
            
        d->G[IJK] -= PORVALNH*fabs(p->W22)*
                    (p->A523*detady + (1.0-p->A523)*detady_n);
        }
        
        if(p->wet[IJp1]==1 && p->wet[IJm1]==0)
        {
        detady = (eta(i,j+1)-eta(i,j))/p->DYP[JP];
        detady_n = (eta_n(i,j+1)-eta_n(i,j))/p->DYP[JP];
            
        d->G[IJK] -= PORVALNH*fabs(p->W22)*
                    (p->A523*detady + (1.0-p->A523)*detady_n);
        }
    }
}

void nhflow_pjm_c::wpgrad(lexer*p, fdm_nhf *d, slice &eta, slice &eta_n)
{
}

double nhflow_pjm_c::limiter(double v1, double v2)
{
    denom = fabs(v1) + fabs(v2);
    
    denom = fabs(denom)>1.0e-10?denom:1.0e10;
    
    val =  (v1*fabs(v2) + fabs(v1)*v2)/denom;

    return val;	
}
