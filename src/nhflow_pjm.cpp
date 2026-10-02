/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#include"nhflow_pjm.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"nhflow_membrane_beta.h"
#include"ghostcell.h"
#include"nhflow_poisson.h"
#include"solver.h"
#include"ioflow.h"
#include"nhflow_poisson.h"
#include"density_f.h"
#include"patchBC_interface.h"
#include"vrans.h"
#include"vrans_definitions.h"

#define WLVL (fabs(WL(i,j))>(p->A544)?WL(i,j):1.0e20)

nhflow_pjm::nhflow_pjm(lexer* p, fdm_nhf *d, ghostcell *pgc, patchBC_interface *ppBC) : teta(1.0)
{
    pBC = ppBC;
    
	pd = new density_f(p);

    ppois = new nhflow_poisson(p);

    gcval_press=540;
    
    gamma=0.5;
    
    p->Darray(P0,p->imax*p->jmax*(p->kmax+2));
}

nhflow_pjm::~nhflow_pjm()
{
}

void nhflow_pjm::start(lexer *p, fdm_nhf *d, solver* psolv, ghostcell* pgc, ioflow *pflow, slice &WL,
                        double *UH, double *VH, double *WH, double alpha)
{
    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<".";
    
    FLOOP
    P0[FIJK] = d->P[FIJK];

    rhs(p,d,pgc,d->U,d->V,d->W,alpha);

    ppois->start(p,d,d->P);

        starttime=pgc->timer();

    psolv->startF(p,pgc,d->P,d->rhsvec,d->M,8);

        endtime=pgc->timer();

	pgc->start7P(p,d->P,gcval_press);

	ucorr(p,d,WL,UH,d->P,alpha);
	vcorr(p,d,WL,VH,d->P,alpha);
	wcorr(p,d,WL,WH,d->P,alpha);
    
    // membranes (X 330): further projection passes. The collocated correction (wide gradient) and the
    // compact matrix differ most where the mobility jumps; repeated passes converge the wide divergence
    // of the corrected velocity to zero (the compact matrix acts as preconditioner).
    for(int it=1; it<d->MPROJ; ++it)
    {
        nhflow_membrane_velupdate(p,d,pgc,WL,UH,VH,WH);
        
        rhs(p,d,pgc,d->U,d->V,d->W,alpha);
        
        FLOOP
        P0[FIJK] = 0.0;
        
        ppois->start(p,d,P0);
        psolv->startF(p,pgc,P0,d->rhsvec,d->M,8);
        pgc->start7P(p,P0,gcval_press);
        
        // the corrections read d->P: apply the increment, then add it to the pressure
        double *Psave = d->P;
        d->P = P0;
        ucorr(p,d,WL,UH,d->P,alpha);
        vcorr(p,d,WL,VH,d->P,alpha);
        wcorr(p,d,WL,WH,d->P,alpha);
        d->P = Psave;
        
        FLOOP
        d->P[FIJK] += P0[FIJK];
        
        pgc->start7P(p,d->P,gcval_press);
    }
    
    // membranes (X 330): face corrections for the Rhie-Chow continuity flux of the next stage
    if(d->MBETA!=nullptr)
    nhflow_membrane_rc_store(p,d,pgc,d->P,alpha*p->dt);

    p->poissoniter=p->solveriter;

	p->poissontime=endtime-starttime;

	if(p->mpirank==0 && p->count%p->P12==0)
	cout<<"piter: "<<p->solveriter<<"  ptime: "<<setprecision(3)<<p->poissontime<<endl;
}

// mesh refinement: start() without the solve (the membrane passes X 330 are not supported there)
double* nhflow_pjm::amr_prepare(lexer *p, fdm_nhf *d, ghostcell *pgc, double alpha)
{
    FLOOP
    P0[FIJK] = d->P[FIJK];
    
    rhs(p,d,pgc,d->U,d->V,d->W,alpha);
    ppois->start(p,d,d->P);
    
    return d->P;
}

void nhflow_pjm::amr_finish(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &WL, double *UH, double *VH, double *WH, double alpha)
{
	pgc->start7P(p,d->P,gcval_press);
    
	ucorr(p,d,WL,UH,d->P,alpha);
	vcorr(p,d,WL,VH,d->P,alpha);
	wcorr(p,d,WL,WH,d->P,alpha);
}

void nhflow_pjm::ucorr(lexer* p, fdm_nhf *d, slice &WL, double *UH, double *P, double alpha)
{
    // membranes (X 330): face mobilities in the horizontal gradient, see nhflow_membrane_beta.h
	LOOP
    WETDRYDEEP
	UH[IJK] -= alpha*p->dt*CPORNH*WL(i,j)*(1.0/p->W1)*
    
                ((d->MBETA!=nullptr ? nhflow_membrane_gradx(p,d,d->P,i,j,k) :
                ((0.5*(d->P[FIp1JKp1]+d->P[FIp1JK])-0.5*(d->P[FIm1JKp1]+d->P[FIm1JK]))/(p->DXP[IP]+p->DXP[IM1])))
                
                + MBETAVAL*0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*((d->P[FIJKp1]-d->P[FIJK])/p->DZN[KP]));
}

void nhflow_pjm::vcorr(lexer* p, fdm_nhf *d, slice &WL, double *VH, double *P, double alpha)
{
    if(p->j_dir==1)
    LOOP
    WETDRYDEEP
    VH[IJK] -= alpha*p->dt*CPORNH*WL(i,j)*(1.0/p->W1)*
    
                ((d->MBETA!=nullptr ? nhflow_membrane_grady(p,d,d->P,i,j,k) :
                ((0.5*(d->P[FIJp1Kp1]+d->P[FIJp1K])-0.5*(d->P[FIJm1Kp1]+d->P[FIJm1K]))/(p->DYP[JP]+p->DYP[JM1])))
                
                + MBETAVAL*0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*((d->P[FIJKp1]-d->P[FIJK])/p->DZN[KP]));
}

void nhflow_pjm::wcorr(lexer* p, fdm_nhf *d, slice &WL, double *WH, double *P, double alpha)
{
    LOOP
    WETDRYDEEP
	WH[IJK] -= MBETAVAL*alpha*p->dt*CPORNH*(1.0/p->W1)*((d->P[FIJKp1]-d->P[FIJK])/(p->DZN[KP]));
}

void nhflow_pjm::rhs(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W, double alpha)
{
    double U1,U2,V1,V2,fac;
    double Um1,Up,Up1;
    double Vp;
    double dz1,dz2;
    double dUdz,dVdz;
    double z,z0,z1,z2;
    double f0,f1,f2;
    double dWdz;
    
    n=0;
    FBASELOOP
    {
	d->rhsvec.V[n]=0.0;
    ++n;
    }

    n=0;
    LOOP
    {
    
        WETDRYDEEP
        {
        fac = p->DZN[KM1]/(p->DZN[KP]+p->DZN[KM1]);    

        U1 = (1.0-fac)*U[Im1JK] + fac*U[Im1JKm1]; 
        U2 = (1.0-fac)*U[Ip1JK] + fac*U[Ip1JKm1]; 
        
        V1 = (1.0-fac)*V[IJm1K] + fac*V[IJm1Km1]; 
        V2 = (1.0-fac)*V[IJp1K] + fac*V[IJp1Km1];     
        
        if(k==0)
        {
        fac = MAX((1.0 - p->A522*fabs(d->Bx(i,j))),0.0)*p->DZN[KM1]/(p->DZN[KP]+p->DZN[KM1]);
        
        U1 = (1.0-fac)*U[Im1JK] + fac*U[Im1JKm1]; 
        U2 = (1.0-fac)*U[Ip1JK] + fac*U[Ip1JKm1]; 
        
        
        fac = MAX((1.0 - p->A522*fabs(d->By(i,j))),0.0)*p->DZN[KM1]/(p->DZN[KP]+p->DZN[KM1]);
        
        V1 = (1.0-fac)*V[IJm1K] + fac*V[IJm1Km1]; 
        V2 = (1.0-fac)*V[IJp1K] + fac*V[IJp1Km1]; 
        }
        

        dUdz = (U[IJK] - U[IJKm1])/p->DZP[KM1];
        dVdz = (V[IJK] - V[IJKm1])/p->DZP[KM1];
        dWdz = p->sigz[IJ]*(W[IJK]-W[IJKm1])/p->DZP[KM1];
        
        d->rhsvec.V[n] =      -  ((U2-U1)/(p->DXP[IP] + p->DXP[IM1])
                                + p->sigx[FIJK]*dUdz
                                
                                + (V2-V1)/(p->DYP[JP] + p->DYP[JM1])
                                + p->sigy[FIJK]*dVdz

                                + dWdz)/(alpha*p->dt);
                                
        }
                            
    ++n;
    }
}

void nhflow_pjm::bedbc(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W,double alpha)
{
}

void nhflow_pjm::vel_setup(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W,double alpha)
{
}

void nhflow_pjm::upgrad(lexer*p, fdm_nhf *d, slice &WL)
{
    LOOP
    WETDRY
    d->F[IJK] += d->eta(i,j)*fabs(p->W22)*
                (d->dfx(i,j) - d->dfx(i-1,j))/(p->DXN[IP]);
}

void nhflow_pjm::vpgrad(lexer*p,fdm_nhf *d, slice &WL)
{
    LOOP
    WETDRY
	d->G[IJK] += d->eta(i,j)*fabs(p->W22)*
                 (d->dfy(i,j) - d->dfy(i,j-1))/(p->DYN[JP]);
}

void nhflow_pjm::wpgrad(lexer*p, fdm_nhf *d, slice &WL)
{
}

