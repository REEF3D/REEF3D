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

#include"nhflow_pjm_hs.h"
#include"lexer.h"
#include"fdm_nhf.h" 
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

#define HX (fabs(d->hx(i,j))>1.0e-20?d->hx(i,j):1.0e20)
#define HXP (fabs(0.5*(d->WL(i,j)+d->WL(i+1,j)))>1.0e-20?0.5*(d->WL(i,j)+d->WL(i+1,j)):1.0e20)
#define HY (fabs(d->hy(i,j))>1.0e-20?d->hy(i,j):1.0e20)
 
nhflow_pjm_hs::nhflow_pjm_hs(lexer* p, fdm_nhf *d, heat *&pheat, concentration *&pconc)
{
    if((p->F80==0||p->A10==5) && p->H10==0 && p->W30==0)
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
}

nhflow_pjm_hs::~nhflow_pjm_hs()
{
}

void nhflow_pjm_hs::start(fdm_nhf *d,lexer*p, poisson* ppois,solver* psolv, ghostcell* pgc, ioflow *pflow, double *U, double *V, double *W, double alpha)
{
}

void nhflow_pjm_hs::ucorr(lexer* p, fdm_nhf *d, double *U, double alpha)
{	
}

void nhflow_pjm_hs::vcorr(lexer* p, fdm_nhf *d, double *U, double alpha)
{	 
}

void nhflow_pjm_hs::wcorr(lexer* p, fdm_nhf *d, double *W, double alpha)
{
}
 
void nhflow_pjm_hs::rhs(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W, double alpha)
{
}
 
void nhflow_pjm_hs::vel_setup(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W, double alpha)
{
}

void nhflow_pjm_hs::upgrad(lexer*p, fdm_nhf *d, slice &eta, slice &eta_n)
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
}

void nhflow_pjm_hs::vpgrad(lexer*p, fdm_nhf *d, slice &eta, slice &eta_n)
{
    if(p->D38==1 && p->A540==1)
    VLOOP
	d->G(i,j,k) -= PORVAL2*fabs(p->W22)*(p->A223*eta(i,j+1) + (1.0-p->A223)*eta_n(i,j+1) - p->A223*eta(i,j) - (1.0-p->A223)*eta_n(i,j))/p->DYP[JP];
    
    if(p->D38==1 && p->A540==2)
    VLOOP
	d->G(i,j,k) -= PORVAL2*fabs(p->W22)*(d->eta(i,j+1) - d->eta(i,j))/p->DYP[JP];
}

void nhflow_pjm_hs::wpgrad(lexer*p, fdm_nhf *d, slice &eta, slice &eta_n)
{
}





