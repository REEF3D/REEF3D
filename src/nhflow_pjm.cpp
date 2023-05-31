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

#include"nhflow_pjm.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"nhflow_poisson.h"
#include"solver.h"
#include"ioflow.h"
#include"nhflow_poisson.h"
#include"density_f.h"
#include"patchBC_interface.h"

nhflow_pjm::nhflow_pjm(lexer* p, fdm_nhf *d, ghostcell *pgc, patchBC_interface *ppBC) : teta(1.0)
{
    pBC = ppBC;
    
	pd = new density_f(p);

    ppois = new nhflow_poisson(p);

    gcval_press=540;
    
    if(p->D33==0)
    solver_id = 8;
    
    if(p->D33==1)
    solver_id = 9;
}

nhflow_pjm::~nhflow_pjm()
{
}

void nhflow_pjm::start(lexer *p, fdm_nhf *d, solver* psolv, ghostcell* pgc, ioflow *pflow,
                        double *U, double *V, double *W, double alpha)
{
    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<".";

    rhs(p,d,pgc,U,V,W,alpha);

    ppois->start(p,d,d->P);

        starttime=pgc->timer();

    psolv->startF(p,pgc,d->P,d->rhsvec,d->M,solver_id);

        endtime=pgc->timer();

	pgc->start7P(p,d->P,gcval_press);
    
    
	ucorr(p,d,U,alpha);
	vcorr(p,d,V,alpha);
	wcorr(p,d,W,alpha);


    p->poissoniter=p->solveriter;

	p->poissontime=endtime-starttime;

	if(p->mpirank==0 && p->count%p->P12==0)
	cout<<"piter: "<<p->solveriter<<"  ptime: "<<setprecision(3)<<p->poissontime<<endl;
}

void nhflow_pjm::ucorr(lexer* p, fdm_nhf *d, double *U, double alpha)
{
	ULOOP
    WETDRY
    if(d->breaking(i,j)==0 && d->breaking(i-1,j)==0 && d->breaking(i+1,j)==0)
    {
	U[IJK] -= alpha*p->dt*CPORNH*PORVALNH*(1.0/p->W1)*
    
                (((0.5*(d->P[FIp1JKp1]+d->P[FIp1JK])-0.5*(d->P[FIJKp1]+d->P[FIJK]))/(p->DXP[IP]))
                
                + 0.5*(p->sigx4[IJK]+p->sigx4[Ip1JK])*((0.5*(d->P[FIJKp1]+d->P[FIp1JKp1])-0.5*(d->P[FIJK]+d->P[FIp1JK]))/p->DZN[KP]));
    }
}

void nhflow_pjm::vcorr(lexer* p, fdm_nhf *d, double *V, double alpha)
{
    VLOOP
    WETDRY
    if(d->breaking(i,j)==0 && d->breaking(i,j-1)==0 && d->breaking(i,j-1)==0)
    {
    V[IJK] -= alpha*p->dt*CPORNH*PORVALNH*(1.0/p->W1)*
    
                (((0.5*(d->P[FIJp1Kp1]+d->P[FIJp1K])-0.5*(d->P[FIJKp1]+d->P[FIJK]))/(p->DYP[JP]))
                
                + 0.5*(p->sigy4[IJK]+p->sigy4[IJp1K])*((0.5*(d->P[FIJKp1]+d->P[FIJp1Kp1])-0.5*(d->P[FIJK]+d->P[FIJp1K]))/p->DZN[KP]));
    }
}

void nhflow_pjm::wcorr(lexer* p, fdm_nhf *d, double *W, double alpha)
{
    WLOOP
    WETDRY
    if(d->breaking(i,j)==0)
    {
    dfdz_plus = (d->P[IJKp2]-d->P[IJKp1])/p->DZN[KP1];
    dfdz_min  = (d->P[IJKp1]-d->P[IJK])/p->DZN[KP];
        
    //detadz = limiter(dfdz_plus,dfdz_min);
    
    detadz = (d->P[IJKp2]-d->P[IJKp1])/(p->DZN[KP1]);
    
	W[IJK] -= alpha*p->dt*CPORNH*PORVALNH*(1.0/p->W1)*detadz*p->sigz[IJ];
    //W[IJK] -= alpha*p->dt*CPORNH*PORVALNH*(1.0/p->W1)*((d->P[FIJK]-d->P[FIJKm1])/(p->DZN[KM1]))*p->sigz[IJ];
    }
}

void nhflow_pjm::rhs(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W, double alpha)
{
    double U1,U2,V1,V2,fac;
    n=0;
    FLOOP
    {
	d->rhsvec.V[n]=0.0;
    ++n;
    }

    n=0;
    LOOP
    {
    fac = p->DZP[KM1]/(p->DZP[KP]+p->DZP[KM1]);
    
    U1 = (1.0-fac)*U[Im1JK] + fac*U[Im1JKm1]; 
    U2 = (1.0-fac)*U[IJK] + fac*U[IJKm1]; 
    
    V1 = (1.0-fac)*V[IJm1K] + fac*V[IJm1Km1]; 
    V2 = (1.0-fac)*V[IJK] + fac*V[IJKm1]; 
    
    dfdz_plus = (W[IJKp1]-W[IJK])/p->DZN[KP];
    dfdz_min  = (W[IJK]-W[IJKm1])/p->DZN[KM1];
        
    //detadz = limiter(dfdz_plus,dfdz_min);
    
    detadz = (W[IJKp1]-W[IJK])/(p->DZN[KP]);
         
    d->rhsvec.V[n] =        -  ((U2-U1)/(p->DXN[IP])
                            + p->sigx[FIJK]*(0.5*(U[IJK]+U[Im1JK]) - 0.5*(U[Im1JKm1]+U[IJKm1]))/p->DZP[KM1]
                            
                            + (V2-V1)/(p->DYN[JP])
                            + p->sigy[FIJK]*(0.5*(V[IJK]+V[IJm1K]) - 0.5*(V[IJm1Km1]+V[IJKm1]))/p->DZP[KM1]

                            + p->sigz[IJ]*detadz)/(alpha*p->dt);
                            
    d->test[IJK] =      -  ((U2-U1)/(p->DXN[IP])
                            + p->sigx[FIJK]*(0.5*(U[IJK]+U[Im1JK]) - 0.5*(U[Im1JKm1]+U[IJKm1]))/p->DZP[KM1]
                            

                            + p->sigz[IJ]*detadz)/(alpha*p->dt);
                            
    ++n;
    }
}

void nhflow_pjm::bedbc(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W,double alpha)
{
}

void nhflow_pjm::vel_setup(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W,double alpha)
{
}

void nhflow_pjm::upgrad(lexer*p, fdm_nhf *d, slice &eta, slice &eta_n)
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

void nhflow_pjm::vpgrad(lexer*p,fdm_nhf *d, slice &eta, slice &eta_n)
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

void nhflow_pjm::wpgrad(lexer*p, fdm_nhf *d, slice &eta, slice &eta_n)
{
}

double nhflow_pjm::limiter(double v1, double v2)
{
    denom = fabs(v1) + fabs(v2);
    
    denom = fabs(denom)>1.0e-10?denom:1.0e10;
    
    val =  (v1*fabs(v2) + fabs(v1)*v2)/denom;

    return val;	
}
