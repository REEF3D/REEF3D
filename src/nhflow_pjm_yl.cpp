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

#define WLVL (fabs(WL(i,j))>(1.0*p->A544)?WL(i,j):1.0e20)

#include"nhflow_pjm_yl.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"nhflow_poisson.h"
#include"solver.h"
#include"ioflow.h"
#include"nhflow_poisson.h"
#include"density_f.h"
#include"patchBC_interface.h"
#include"vrans.h"

nhflow_pjm_yl::nhflow_pjm_yl(lexer* p, fdm_nhf *d, ghostcell *pgc, patchBC_interface *ppBC) : teta(1.0)
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
    
    gamma=0.5;
}

nhflow_pjm_yl::~nhflow_pjm_yl()
{
}

void nhflow_pjm_yl::start(lexer *p, fdm_nhf *d, solver* psolv, ghostcell* pgc, ioflow *pflow, slice &WL,
                        double *UH, double *VH, double *WH, double alpha)
{
    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<".";
    
        starttime=pgc->timer();

    press_integral(p,d,pgc,d->U,d->V,d->W,alpha);
    
    press_corr(p,d,WL,d->P,PCORR,alpha);
    
    pgc->start7P(p,d->P,gcval_press);
    pgc->start7P(p,PCORR,gcval_press);
    
	ucorr(p,d,WL,UH,PCORR,alpha);
	vcorr(p,d,WL,VH,PCORR,alpha);
	wcorr(p,d,WL,WH,PCORR,alpha);

    p->poissoniter=p->solveriter;
    
        endtime=pgc->timer();
    
	p->ptime=endtime-starttime;
    p->poissontime+=p->ptime;
    

	if(p->mpirank==0 && p->count%p->P12==0)
	cout<<"piter: "<<p->solveriter<<"  ptime: "<<setprecision(3)<<p->ptime<<endl;
}

void nhflow_pjm_yl::press_corr(lexer* p, fdm_nhf *d, slice &WL, double *P, double *PCORR, double alpha)
{
	FLOOP
    WETDRYDEEP
    PCORR[FIJK] = P[FIJK] - PCORR[FIJK];
}

void nhflow_pjm_yl::press_integral(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W, double alpha)
{
    FLOOP
    WETDRYDEEP
    PCORR[FIJK] = d->P[FIJK];
    
    LOOP
    d->test[IJK] = 0.0;
    
    k=p->knoz;
    SLICELOOP4
    {
    d->P[FIJK]=0.0;
    ++n;
    }
    
    pgc->start7P(p,PCORR,gcval_press);
    
    YLLOOP
    WETDRYDEEP
    {
    d->P[FIJK] =      d->P[FIJKp1] 
    
                    + p->WL[IJ]*p->DZN[KP]*d->DWDT[IJK]  
                    
                    + p->WL[IJ]*(W[IJK]-W[IJKm1])*p->sigt[FIJK]
                    
                    + p->DZN[KP]*d->FSW[IJK];
                    
                    
    d->test[IJK] =    d->test[IJKp1]
                    
                    + p->WL[IJ]*p->DZN[KP]*d->DWDT[IJK]  ;
                    
    }
    
    pgc->start7P(p,PCORR,gcval_press);
}

void nhflow_pjm_yl::bedbc(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W,double alpha)
{
}

void nhflow_pjm_yl::vel_setup(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W,double alpha)
{
}

void nhflow_pjm_yl::ucorr(lexer* p, fdm_nhf *d, slice &WL, double *UH, double *PCORR, double alpha)
{
	LOOP
    WETDRYDEEP
    {
	UH[IJK] -= alpha*p->dt*WL(i,j)*(1.0/p->W1)*
    
                ((0.5*(PCORR[FIp1JKp1]+PCORR[FIp1JK])-0.5*(PCORR[FIm1JKp1]+PCORR[FIm1JK]))/(p->DXP[IP]+p->DXP[IM1])
                
                + 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*((PCORR[FIJKp1]-PCORR[FIJK])/p->DZN[KP]));
    }
}

void nhflow_pjm_yl::vcorr(lexer* p, fdm_nhf *d, slice &WL, double *VH, double *PCORR, double alpha)
{
    if(p->j_dir==1)
    LOOP
    WETDRYDEEP
    {
    VH[IJK] -= alpha*p->dt*WL(i,j)*(1.0/p->W1)*
    
                ((0.5*(PCORR[FIJp1Kp1]+PCORR[FIJp1K])-0.5*(PCORR[FIJm1Kp1]+PCORR[FIJm1K]))/(p->DYP[JP]+p->DYP[JM1])
                
                + 0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*((PCORR[FIJKp1]-PCORR[FIJK])/p->DZN[KP]));
    }
}

void nhflow_pjm_yl::wcorr(lexer* p, fdm_nhf *d, slice &WL, double *WH, double *PCORR, double alpha)
{
    /*LOOP
    WETDRYDEEP
	WH[IJK] -= alpha*p->dt*(1.0/p->W1)*((PCORR[FIJKp1]-PCORR[FIJK])/(p->DZN[KP]));*/
}

void nhflow_pjm_yl::upgrad(lexer*p, fdm_nhf *d, slice &WL)
{
    LOOP
    WETDRY
    d->F[IJK] += d->eta(i,j)*fabs(p->W22)*
                (d->dfx(i,j) - d->dfx(i-1,j))/(p->DXN[IP]);
                
    LOOP
    WETDRYDEEP
    {
    dPdx = (0.5*(d->P[FIp1JKp1]+d->P[FIp1JK])-0.5*(d->P[FIm1JKp1]+d->P[FIm1JK]))/(p->DXP[IP]+p->DXP[IM1]);
    
    d->F[IJK] -= (1.0/p->W1)*WL(i,j)*
                (dPdx
                + 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*((d->P[FIJKp1]-d->P[FIJK])/p->DZN[KP]));
    }
}

void nhflow_pjm_yl::vpgrad(lexer*p, fdm_nhf *d, slice &WL)
{
    
    if(p->j_dir==1)
    LOOP
    WETDRY
	d->G[IJK] += d->eta(i,j)*fabs(p->W22)*
                 (d->dfy(i,j) - d->dfy(i,j-1))/(p->DYN[JP]);
       
    if(p->j_dir==1) 
    LOOP
    WETDRYDEEP
    {
    dPdy = (0.5*(d->P[FIJp1Kp1]+d->P[FIJp1K])-0.5*(d->P[FIJm1Kp1]+d->P[FIJm1K]))/(p->DYP[JP]+p->DYP[JM1]);
    
	d->G[IJK] -= (1.0/p->W1)*WL(i,j)*
                (dPdy
                + 0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*((d->P[FIJKp1]-d->P[FIJK])/p->DZN[KP]));
    }
}

void nhflow_pjm_yl::wpgrad(lexer*p, fdm_nhf *d, slice &WL)
{
    LOOP
    WETDRYDEEP
    d->H[IJK] -= (1.0/p->W1)*((d->P[FIJKp1]-d->P[FIJK])/(p->DZN[KP]));
}


