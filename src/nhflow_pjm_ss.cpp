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

#include"nhflow_pjm_ss.h"
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
#include"hypre_struct.h"
#include"hypre_sstruct_fnpf.h"
 
nhflow_pjm_ss::nhflow_pjm_ss(lexer* p, fdm_nhf *d, ghostcell *pgc)
{
	pd = new density_f(p);

    gcval_press=540;  

	gcval_u=10;
	gcval_v=11;
	gcval_w=12;

    vecsize=p->knox*p->knoy*(p->knoz+1); 
    
    p->Darray(M,vecsize*15);
    p->Darray(x,vecsize);
    p->Darray(rhs,vecsize);
    
    teta=0.5;
}

nhflow_pjm_ss::~nhflow_pjm_ss()
{
}

void nhflow_pjm_ss::start(lexer*p, fdm_nhf *d, solver* psolv, ghostcell* pgc, ioflow *pflow, double *U, double *V, double *W, double alpha)
{    
    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<".";
				
    rhscalc(p,d,pgc,U,V,W,alpha);
    
    if(p->j_dir==0)
    poisson2D(p,d,d->press);
    
    if(p->j_dir==1)
    poisson3D(p,d,d->press);
    
    fillvec(p,d,pgc);
	
        starttime=pgc->timer();

    psolv->startM(p,pgc,x,rhs,M,5);

        endtime=pgc->timer();
        
    fillvec_back(p,d,pgc);
    
	pgc->start4(p,d->press,gcval_press);
    
	ucorr(p,d,U,alpha);
	vcorr(p,d,V,alpha);
	wcorr(p,d,W,alpha);

    p->poissoniter=p->solveriter;

	p->poissontime=endtime-starttime;

	if(p->mpirank==0 && p->count%p->P12==0)
	cout<<"piter: "<<p->solveriter<<"  ptime: "<<setprecision(3)<<p->poissontime<<endl;
    
    //pgc->start4(p,d->test,1);
}

void nhflow_pjm_ss::ucorr(lexer* p, fdm_nhf *d, double *U, double alpha)
{	
	if(p->D37==1)
	LOOP
	U[IJK] -= alpha*p->dt*CPORNH*PORVALNH*(1.0/p->W1)*((d->press(i+1,j,k)-d->press(i,j,k))/p->DXP[IP]
                + 0.25*(p->sigx[FIJK]+p->sigx[FIJKp1]+p->sigx[FIp1JK]+p->sigx[FIp1JKp1])*(0.5*(d->press(i,j,k+1)+d->press(i+1,j,k+1))-0.5*(d->press(i,j,k-1)+d->press(i+1,j,k-1)))/(p->DZP[KP]+p->DZP[KP1]));
    
    if(p->D37==2)
    LOOP
    {       
     check=0;
    
        if(p->flag1[IJKp1]<0)
        check=1;        
    
    if(check==1)
    U[IJK] -= alpha*p->dt*CPORNH*PORVALNH*(1.0/p->W1)*((d->press(i+1,j,k)-d->press(i,j,k))/p->DXP[IP]
                + 0.25*(p->sigx[FIJK]+p->sigx[FIJKp1]+p->sigx[FIp1JK]+p->sigx[FIp1JKp1])*(0.5*((1.0 - 1.0/teta)*(d->press(i,j,k)+d->press(i+1,j,k)))-0.5*(d->press(i,j,k-1)+d->press(i+1,j,k-1)))/(p->DZP[KP]+p->DZP[KP1]));
    
    if(check==0)
    U[IJK] -= alpha*p->dt*CPORNH*PORVALNH*(1.0/p->W1)*((d->press(i+1,j,k)-d->press(i,j,k))/p->DXP[IP]
                + 0.25*(p->sigx[FIJK]+p->sigx[FIJKp1]+p->sigx[FIp1JK]+p->sigx[FIp1JKp1])*(0.5*(d->press(i,j,k+1)+d->press(i+1,j,k+1))-0.5*(d->press(i,j,k-1)+d->press(i+1,j,k-1)))/(p->DZP[KP]+p->DZP[KP1]));
    }
}

void nhflow_pjm_ss::vcorr(lexer* p, fdm_nhf *d, double *V, double alpha)
{	 
    if(p->D37==1)
    LOOP
    V[IJK] -= alpha*p->dt*CPORNH*PORVALNH*(1.0/p->W1)*((d->press(i,j+1,k)-d->press(i,j,k))/p->DYP[JP] 
                + 0.25*(p->sigy[FIJK]+p->sigy[FIJKp1]+p->sigy[FIJp1K]+p->sigy[FIJp1Kp1])*(0.5*(d->press(i,j,k+1)+d->press(i,j+1,k+1))-0.5*(d->press(i,j,k-1)+d->press(i,j+1,k-1)))/(p->DZP[KP]+p->DZP[KP1]));
                
                
    if(p->D37==2)
    LOOP
    {       
     check=0;
    
        if(p->flag2[IJKp1]<0)
        check=1;        
    
    if(check==1)
    V[IJK] -= alpha*p->dt*CPORNH*PORVALNH*(1.0/p->W1)*((d->press(i,j+1,k)-d->press(i,j,k))/p->DYP[JP] 
                + 0.25*(p->sigy[FIJK]+p->sigy[FIJKp1]+p->sigy[FIJp1K]+p->sigy[FIJp1Kp1])*(0.5*(1.0 - 1.0/teta)*(d->press(i,j,k)+d->press(i,j+1,k))-0.5*(d->press(i,j,k-1)+d->press(i,j+1,k-1)))/(p->DZP[KP]+p->DZP[KP1]));
    
    if(check==0)
    V[IJK] -= alpha*p->dt*CPORNH*PORVALNH*(1.0/p->W1)*((d->press(i,j+1,k)-d->press(i,j,k))/p->DYP[JP] 
                + 0.25*(p->sigy[FIJK]+p->sigy[FIJKp1]+p->sigy[FIJp1K]+p->sigy[FIJp1Kp1])*(0.5*(d->press(i,j,k+1)+d->press(i,j+1,k+1))-0.5*(d->press(i,j,k-1)+d->press(i,j+1,k-1)))/(p->DZP[KP]+p->DZP[KP1]));
    }
}

void nhflow_pjm_ss::wcorr(lexer* p, fdm_nhf *d, double *W, double alpha)
{
    if(p->D37==1)
    LOOP 	
	W[IJK] -= alpha*p->dt*CPORNH*PORVALNH*((d->press(i,j,k+1)-d->press(i,j,k))/(p->DZP[KP]*p->W1))*p->sigz[IJ];
    
    
    if(p->D37==2)
	LOOP
    {
    check=0;
    
        if(p->flag3[IJKp1]<0)
        check=1;

    if(check==1)    
    W[IJK] -= alpha*p->dt*CPORNH*PORVALNH*(((1.0 - 1.0/teta)*d->press(i,j,k)-d->press(i,j,k))/(p->DZP[KP]*p->W1))*p->sigz[IJ];
    
    if(check==0)
    W[IJK] -= alpha*p->dt*CPORNH*PORVALNH*((d->press(i,j,k+1)-d->press(i,j,k))/(p->DZP[KP]*p->W1))*p->sigz[IJ];
    }
}
 
void nhflow_pjm_ss::rhscalc(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W, double alpha)
{
    NLOOP4
	rhs[n]=0.0;
	
    pip=p->Y50;

    n=0;
    KJILOOP
    {
    PCHECK
    rhs[n] =            -  ((U[IJK]-U[Im1JK])/p->DXN[IP]
                            + 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*(0.5*(U[IJKp1]+U[Im1JKp1])-0.5*(U[IJKm1]+U[Im1JKm1]))/(p->DZP[KP]+p->DZP[KP1])
                            
                            + (V[IJK]-U[IJm1K])/p->DYN[JP] 
                            + 0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*(0.5*(V[IJKp1]+V[IJm1Kp1])-0.5*(V[IJKm1]+V[IJm1Km1]))/(p->DZP[KP]+p->DZP[KP1])
                           
                            + p->sigz[IJ]*(W[IJK]-W[IJKm1])/p->DZN[KP])/(alpha*p->dt);
                           
    SCHECK
    rhs[n] = 0.0;
                                        
    ++n;
    }
    pip=0;
}
 
void nhflow_pjm_ss::vel_setup(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W, double alpha)
{
	pgc->start7V(p,U,d->bc,gcval_u);    pgc->start7V(p,V,d->bc,gcval_v);    pgc->start7V(p,W,d->bc,gcval_w);
}

void nhflow_pjm_ss::upgrad(lexer*p,fdm_nhf *d, slice &eta, slice &eta_n)
{
    if(p->D38==1 && p->A540==1)
    LOOP
	d->F[IJK] -= PORVALNH*fabs(p->W22)*(p->A223*eta(i+1,j) + (1.0-p->A223)*eta_n(i+1,j) - p->A223*eta(i,j) - (1.0-p->A223)*eta_n(i,j))/p->DXP[IP];
    
    if(p->D38==1 && p->A540==2)
    LOOP
	d->F[IJK] -= PORVALNH*fabs(p->W22)*(d->eta(i+1,j) - d->eta(i,j))/p->DXP[IP];
    
   
    if(p->D38==2 && p->A540==1)
    LOOP
	d->F[IJK] -= PORVALNH*fabs(p->W22)*(1.0/HX)*
    
                    (0.5*(pow(eta(i+1,j),2.0) - pow(eta(i,j),2.0))/p->DXP[IP]
                    
                    + ((p->A223*eta(i+1,j) + (1.0-p->A223)*eta_n(i+1,j))*d->depth(i+1,j) - (p->A223*eta(i,j) + (1.0-p->A223)*eta_n(i,j))*d->depth(i,j))/p->DXP[IP]
                    
                    - 0.5*((p->A223*eta(i,j) + (1.0-p->A223)*eta_n(i,j)) + (p->A223*eta(i+1,j) + (1.0-p->A223)*eta_n(i+1,j)))*(d->depth(i+1,j)-d->depth(i,j))/p->DXP[IP]);
    
    if(p->D38==2 && p->A540==2)
    LOOP
	d->F[IJK] -= PORVALNH*fabs(p->W22)*(1.0/HX)*
    
                    (0.5*(pow(d->eta(i+1,j),2.0) - pow(d->eta(i,j),2.0))/p->DXP[IP]
                    
                    + (d->eta(i+1,j)*d->depth(i+1,j) - d->eta(i,j)*d->depth(i,j))/p->DXP[IP]
                    
                    - 0.5*(d->eta(i,j) + d->eta(i+1,j))*(d->depth(i+1,j)-d->depth(i,j))/p->DXP[IP]);
    
    // fx = 1/2 g (eta^2 - 2* eta *z_b)
    // Sx = -g * eta * eta * Bx
}

void nhflow_pjm_ss::vpgrad(lexer*p,fdm_nhf *d, slice &eta, slice &eta_n)
{
    if(p->D38==1 && p->A540==1)
    LOOP
	d->G[IJK] -= PORVALNH*fabs(p->W22)*(p->A223*eta(i,j+1) + (1.0-p->A223)*eta_n(i,j+1) - p->A223*eta(i,j) - (1.0-p->A223)*eta_n(i,j))/p->DYP[JP];
    
    if(p->D38==1 && p->A540==2)
    LOOP
	d->G[IJK] -= PORVALNH*fabs(p->W22)*(d->eta(i,j+1) - d->eta(i,j))/p->DYP[JP];
}

void nhflow_pjm_ss::wpgrad(lexer*p,fdm_nhf *d, slice &eta, slice &eta_n)
{
}

void nhflow_pjm_ss::fillvec(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    n=0;
    KJILOOP
    {
    PCHECK
    x[n] = d->press(i,j,k);
    
    SCHECK
    x[n] = 0.0;
    
    ++n;
    }
}

void nhflow_pjm_ss::fillvec_back(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    n=0;
    KJILOOP
    {
    PCHECK
    d->press(i,j,k) = x[n];
    
    ++n;
    }
}






