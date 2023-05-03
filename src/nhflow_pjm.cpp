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

#include"nhflow_pjm.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"nhflow_poisson.h"
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

void nhflow_pjm::start(lexer *p, fdm_nhf *d, solver* psolv, ghostcell* pgc, ioflow *pflow,
                        double *U, double *V, double *W, double alpha)
{
    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<".";

    rhs(p,d,pgc,U,V,W,alpha);

    bedbc(p,d,pgc,U,V,W,alpha);

    ppois->start(p,d,d->P);

        starttime=pgc->timer();

    psolv->startF(p,pgc,d->P,d->rhsvec,d->M,8);

        endtime=pgc->timer();

	pgc->start4V(p,d->P,d->bc,gcval_press);

	ucorr(p,d,U,alpha);
	vcorr(p,d,V,alpha);
	wcorr(p,d,W,alpha);

    p->poissoniter=p->solveriter;

	p->poissontime=endtime-starttime;

	if(p->mpirank==0 && p->count%p->P12==0)
	cout<<"piter: "<<p->solveriter<<"  ptime: "<<setprecision(3)<<p->poissontime<<endl;

    //pgc->start4(p,d->test,1);
}

void nhflow_pjm::ucorr(lexer* p, fdm_nhf *d, double *U, double alpha)
{
	LOOP
	U[IJK] -= alpha*p->dt*CPORNH*PORVALNH*(1.0/p->W1)*((d->P[FIp1JK]-d->P[FIJK])/p->DXP[IP]
                + 0.25*(p->sigx[FIJK]+p->sigx[FIJKp1]+p->sigx[FIp1JK]+p->sigx[FIp1JKp1])*(0.5*(d->P[FIJKp1]+d->P[FIp1JKp1])-0.5*(d->P[FIJKm1]+d->P[FIp1JKm1]))/(p->DZP[KP]+p->DZP[KP1]));
}

void nhflow_pjm::vcorr(lexer* p, fdm_nhf *d, double *V,double alpha)
{
    LOOP
    V[IJK] -= alpha*p->dt*CPORNH*PORVALNH*(1.0/p->W1)*((d->P[FIJp1K]-d->P[FIJK])/p->DYP[JP]
                + 0.25*(p->sigy[FIJK]+p->sigy[FIJKp1]+p->sigy[FIJp1K]+p->sigy[FIJp1Kp1])*(0.5*(d->P[FIJKp1]+d->P[FIJp1Kp1])-0.5*(d->P[FIJKm1]+d->P[FIJp1Km1]))/(p->DZP[KP]+p->DZP[KP1]));
}

void nhflow_pjm::wcorr(lexer* p, fdm_nhf *d, double *W, double alpha)
{
    LOOP
	W[IJK] -= alpha*p->dt*CPORNH*PORVALNH*((d->P[FIJKp1]-d->P[FIJK])/(p->DZP[KP]*p->W1))*p->sigz[IJ];
}

void nhflow_pjm::rhs(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W, double alpha)
{
    NLOOP4
	d->rhsvec.V[n]=0.0;

    pip=p->Y50;

    n=0;
    LOOP
    {
    d->rhsvec.V[n] =      -  ((U[IJK]-U[Im1JK])/p->DXN[IP]
                            + 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*(0.5*(U[IJKp1]+U[Im1JKp1])-0.5*(U[IJKm1]+U[Im1JKm1]))/(p->DZP[KP]+p->DZP[KP1])

                            + (V[IJK]-V[IJm1K])/p->DYN[JP]
                            + 0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*(0.5*(V[IJKp1]+V[IJm1Kp1])-0.5*(V[IJKm1]+V[IJm1Km1]))/(p->DZP[KP]+p->DZP[KP1])

                            + p->sigz[IJ]*(W[IJK]-W[IJm1K])/p->DZN[KP])/(alpha*p->dt);

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
	d->F[IJK] -= PORVALNH*fabs(p->W22)*(p->A223*eta(i+1,j) + (1.0-p->A223)*eta_n(i+1,j) - p->A223*eta(i,j) - (1.0-p->A223)*eta_n(i,j))/p->DXP[IP];

    if(p->D38==1 && p->A540==2)
    ULOOP
	d->F[IJK] -= PORVALNH*fabs(p->W22)*(d->eta(i+1,j) - d->eta(i,j))/p->DXP[IP];


    if(p->D38==2 && p->A540==1)
    ULOOP
	d->F[IJK] -= PORVALNH*fabs(p->W22)*(1.0/HX)*

                    (0.5*(pow(eta(i+1,j),2.0) - pow(eta(i,j),2.0))/p->DXP[IP]

                    + ((p->A223*eta(i+1,j) + (1.0-p->A223)*eta_n(i+1,j))*d->depth(i+1,j) - (p->A223*eta(i,j) + (1.0-p->A223)*eta_n(i,j))*d->depth(i,j))/p->DXP[IP]

                    - 0.5*((p->A223*eta(i,j) + (1.0-p->A223)*eta_n(i,j)) + (p->A223*eta(i+1,j) + (1.0-p->A223)*eta_n(i+1,j)))*(d->depth(i+1,j)-d->depth(i,j))/p->DXP[IP]);

    if(p->D38==2 && p->A540==2)
    ULOOP
	d->F[IJK] -= PORVALNH*fabs(p->W22)*(1.0/HX)*

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
	d->G[IJK] -= PORVALNH*fabs(p->W22)*(p->A223*eta(i,j+1) + (1.0-p->A223)*eta_n(i,j+1) - p->A223*eta(i,j) - (1.0-p->A223)*eta_n(i,j))/p->DYP[JP];

    if(p->D38==1 && p->A540==2)
    VLOOP
	d->G[IJK] -= PORVALNH*fabs(p->W22)*(d->eta(i,j+1) - d->eta(i,j))/p->DYP[JP];
}

void nhflow_pjm::wpgrad(lexer*p, fdm_nhf *d, slice &eta, slice &eta_n)
{
}
