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

#define WLVL (fabs(WL(i,j))>(1.0*p->A544)?WL(i,j):1.0e20)

#include"nhflow_pjm_corr.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"nhflow_poisson.h"
#include"solver.h"
#include"ioflow.h"
#include"nhflow_poisson.h"
#include"density_f.h"
#include"patchBC_interface.h"

nhflow_pjm_corr::nhflow_pjm_corr(lexer* p, fdm_nhf *d, ghostcell *pgc, patchBC_interface *ppBC) : teta(1.0)
{
    pBC = ppBC;
    
	pd = new density_f(p);

    ppois = new nhflow_poisson(p);
    
    p->Darray(PCORR,p->imax*p->jmax*(p->kmax+2));
    
    gcval_press=540;
    
    if(p->D33==0)
    solver_id = 8;
    
    if(p->D33==1)
    solver_id = 9;
}

nhflow_pjm_corr::~nhflow_pjm_corr()
{
}

void nhflow_pjm_corr::start(lexer *p, fdm_nhf *d, solver* psolv, ghostcell* pgc, ioflow *pflow, slice &WL,
                        double *UH, double *VH, double *WH, double alpha)
{
    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<".";

    rhs(p,d,pgc,d->U,d->V,d->W,alpha);

    ppois->start(p,d,PCORR);

        starttime=pgc->timer();

    psolv->startF(p,pgc,PCORR,d->rhsvec,d->M,solver_id);

        endtime=pgc->timer();
    
    presscorr(p,d,WL,d->P,PCORR,alpha);
    
    pgc->start7P(p,d->P,gcval_press);
    pgc->start7P(p,PCORR,gcval_press);
    
	ucorr(p,d,WL,UH,PCORR,alpha);
	vcorr(p,d,WL,VH,PCORR,alpha);
	wcorr(p,d,WL,WH,PCORR,alpha);

    p->poissoniter=p->solveriter;

	p->poissontime=endtime-starttime;

	if(p->mpirank==0 && p->count%p->P12==0)
	cout<<"piter: "<<p->solveriter<<"  ptime: "<<setprecision(3)<<p->poissontime<<endl;
}

void nhflow_pjm_corr::presscorr(lexer* p, fdm_nhf *d, slice &WL, double *P, double *PCORR, double alpha)
{
	FLOOP
    WETDRYDEEP
    if(d->breaking(i,j)==0)
    P[FIJK] += PCORR[FIJK];
}

void nhflow_pjm_corr::ucorr(lexer* p, fdm_nhf *d, slice &WL, double *UH, double *PCORR, double alpha)
{
	LOOP
    WETDRYDEEP
    if(d->breaking(i,j)==0 && d->breaking(i-1,j)==0 && d->breaking(i+1,j)==0)
	UH[IJK] -= alpha*p->dt*CPORNH*PORVALNH*WL(i,j)*(1.0/p->W1)*
    
                (((0.5*(PCORR[FIp1JKp1]+PCORR[FIp1JK])-0.5*(PCORR[FIm1JKp1]+PCORR[FIm1JK]))/(p->DXP[IP]+p->DXP[IM1]))
                
                + 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*((PCORR[FIJKp1]-PCORR[FIJK])/p->DZN[KP]));
}

void nhflow_pjm_corr::vcorr(lexer* p, fdm_nhf *d, slice &WL, double *VH, double *PCORR, double alpha)
{
    if(p->j_dir==1)
    LOOP
    WETDRYDEEP
    if(d->breaking(i,j)==0 && d->breaking(i,j-1)==0 && d->breaking(i,j+1)==0)
    VH[IJK] -= alpha*p->dt*CPORNH*PORVALNH*WL(i,j)*(1.0/p->W1)*
    
                (((0.5*(PCORR[FIJp1Kp1]+PCORR[FIJp1K])-0.5*(PCORR[FIJm1Kp1]+PCORR[FIJm1K]))/(p->DYP[JP]+p->DYP[JM1]))
                
                + 0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*((PCORR[FIJKp1]-PCORR[FIJK])/p->DZN[KP]));
}

void nhflow_pjm_corr::wcorr(lexer* p, fdm_nhf *d, slice &WL, double *WH, double *PCORR, double alpha)
{
    LOOP
    WETDRYDEEP
    if(d->breaking(i,j)==0)
	WH[IJK] -= alpha*p->dt*CPORNH*PORVALNH*(1.0/p->W1)*((PCORR[FIJKp1]-PCORR[FIJK])/(p->DZN[KP]));
}

/*
void nhflow_pjm_corr::rhs(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W, double alpha)
{
    double U1,U2,V1,V2,fac;

    n=0;
    FLOOP
    {
	d->rhsvec.V[n]=0.0;
    PCORR[FIJK]=0.0;
    ++n;
    }

    n=0;
    LOOP
    {
    fac = p->DZN[KM1]/(p->DZN[KP]+p->DZN[KM1]);    
    //fac = p->DZN[KP]/(p->DZN[KP]+p->DZN[KM1]);

    U1 = (1.0-fac)*U[Im1JK] + fac*U[Im1JKm1]; 
    U2 = (1.0-fac)*U[Ip1JK] + fac*U[Ip1JKm1]; 
    
    V1 = (1.0-fac)*V[IJm1K] + fac*V[IJm1Km1]; 
    V2 = (1.0-fac)*V[IJp1K] + fac*V[IJp1Km1];     
         
    d->rhsvec.V[n] =      -  ((U2-U1)/(p->DXP[IP] + p->DXP[IM1])
                            + p->sigx[FIJK]*(U[IJK]-U[IJKm1])/p->DZP[KM1]
                            
                            + (V2-V1)/(p->DYP[JP] + p->DYP[JM1])
                            + p->sigy[FIJK]*(V[IJK]-V[IJKm1])/p->DZP[KM1]

                            + p->sigz[IJ]*(W[IJK]-W[IJKm1])/p->DZP[KM1])/(alpha*p->dt);
                            
    ++n;
    }
}*/


void nhflow_pjm_corr::rhs(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W, double alpha)
{
    double U1,U2,V1,V2,fac;
    double Um1,Up,Up1;
    double Vp;
    double dz1,dz2;
    double dUdz,dVdz;
    double z,z0,z1,z2;
    double f0,f1,f2;
    
    n=0;
    FLOOP
    {
	d->rhsvec.V[n]=0.0;
    PCORR[FIJK]=0.0;
    ++n;
    }

    n=0;
    LOOP
    {
    fac = p->DZN[KM1]/(p->DZN[KP]+p->DZN[KM1]);    

    U1 = (1.0-fac)*U[Im1JK] + fac*U[Im1JKm1]; 
    U2 = (1.0-fac)*U[Ip1JK] + fac*U[Ip1JKm1]; 
    
    V1 = (1.0-fac)*V[IJm1K] + fac*V[IJm1Km1]; 
    V2 = (1.0-fac)*V[IJp1K] + fac*V[IJp1Km1];   

/*
    if(k==0)
    {
    U1 = U[Im1JK];
    U2 = U[Ip1JK];
    }*/

 
    // dz
    z0 = p->ZP[KM2];
    z1 = p->ZP[KM1];
    z2 = p->ZP[KP];
    z  = p->ZP[KP] - p->DZN[KP];
    
    f0 = U[IJKm2];
    f1 = U[IJKm1];
    f2 = U[IJK];
    
    Up = f0*(z-z1)*(z-z2)/((z0-z1)*(z0-z2)) + f1*(z-z0)*(z-z2)/((z1-z0)*(z1-z2)) + f2*(z-z0)*(z-z1)/((z2-z0)*(z2-z1));
    
    
    f0 = V[IJKm2];
    f1 = V[IJKm1];
    f2 = V[IJK];
    
    Vp = f0*(z-z1)*(z-z2)/((z0-z1)*(z0-z2)) + f1*(z-z0)*(z-z2)/((z1-z0)*(z1-z2)) + f2*(z-z0)*(z-z1)/((z2-z0)*(z2-z1));
    
    /*
    z0 = p->ZP[KM1];
    z1 = p->ZP[KP];
    z  = p->ZP[KP] - p->DZN[KP];
    
    f0 = U[IJKm1];
    f1 = U[IJK];

    
    Up = f0*(z-z1)/(z0-z1) + f1*(z-z0)/(z1-z0);*/
    
    dUdz = (U[IJK] - Up)/p->DZN[KP];
    
    dVdz = (V[IJK] - Vp)/p->DZN[KP];
    
    //dUdz = (U[IJK] - U[IJKm1])/p->DZN[KP];
     
    d->rhsvec.V[n] =      -  ((U2-U1)/(p->DXP[IP] + p->DXP[IM1])
                            + p->sigx[FIJK]*dUdz
                            
                            + (V2-V1)/(p->DYP[JP] + p->DYP[JM1])
                            + p->sigy[FIJK]*dVdz

                            + p->sigz[IJ]*(W[IJK]-W[IJKm1])/p->DZP[KM1])/(alpha*p->dt);
                            
    ++n;
    }
}

void nhflow_pjm_corr::bedbc(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W,double alpha)
{
}

void nhflow_pjm_corr::vel_setup(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W,double alpha)
{
}

void nhflow_pjm_corr::upgrad(lexer*p, fdm_nhf *d, slice &WL)
{
    LOOP
    //WETDRY
    d->F[IJK] += PORVALNH*0.5*(d->ETAs(i,j)+d->ETAn(i-1,j))*fabs(p->W22)*
                (d->dfx(i,j) - d->dfx(i-1,j))/(p->DXN[IP]);
                
    LOOP
    WETDRYDEEP
    if(d->breaking(i,j)==0 && d->breaking(i-1,j)==0 && d->breaking(i+1,j)==0)
    d->F[IJK] -= PORVALNH*(1.0/p->W1)*WL(i,j)*
                (((0.5*(d->P[FIp1JKp1]+d->P[FIp1JK])-0.5*(d->P[FIm1JKp1]+d->P[FIm1JK]))/(p->DXP[IP]+p->DXP[IM1]))
                + 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*((d->P[FIJKp1]-d->P[FIJK])/p->DZN[KP]));
}

void nhflow_pjm_corr::vpgrad(lexer*p, fdm_nhf *d, slice &WL)
{
    LOOP
    //WETDRY
	d->G[IJK] += PORVALNH*0.5*(d->ETAe(i,j)+d->ETAw(i,j-1))*fabs(p->W22)*
                 (d->dfy(i,j) - d->dfy(i,j-1))/(p->DYN[JP]);
                
    LOOP
    WETDRYDEEP
    if(d->breaking(i,j)==0 && d->breaking(i,j-1)==0 && d->breaking(i,j+1)==0)
	d->G[IJK] -= PORVALNH*(1.0/p->W1)*WL(i,j)*
                (((0.5*(d->P[FIJp1Kp1]+d->P[FIJp1K])-0.5*(d->P[FIJm1Kp1]+d->P[FIJm1K]))/(p->DYP[JP]+p->DYP[JM1]))
                + 0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*((d->P[FIJKp1]-d->P[FIJK])/p->DZN[KP]));
}

void nhflow_pjm_corr::wpgrad(lexer*p, fdm_nhf *d, slice &WL)
{
    LOOP
    WETDRYDEEP
    if(d->breaking(i,j)==0 && d->breaking(i-1,j)==0 && d->breaking(i+1,j)==0 && d->breaking(i,j-1)==0 && d->breaking(i,j+1)==0)
    d->H[IJK] -= PORVALNH*(1.0/p->W1)*((d->P[FIJKp1]-d->P[FIJK])/(p->DZN[KP]));
}

