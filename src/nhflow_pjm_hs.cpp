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

#include"nhflow_pjm_hs.h"
#include"lexer.h"
#include"fdm_nhf.h" 
#include"ghostcell.h"
#include"nhflow_poisson.h"
#include"solver.h"
#include"ioflow.h"
#include"heat.h"
#include"concentration.h"
#include"density_f.h"
#include"patchBC_interface.h"

#define HX (fabs(d->hx(i,j))>1.0e-20?d->hx(i,j):1.0e20)
#define HXP (fabs(0.5*(d->WL(i,j)+d->WL(i+1,j)))>1.0e-20?0.5*(d->WL(i,j)+d->WL(i+1,j)):1.0e20)
#define HY (fabs(d->hy(i,j))>1.0e-20?d->hy(i,j):1.0e20)
#define WLVL (fabs(d->WL(i,j))>1.0e-20?d->WL(i,j):1.0e20)
 
nhflow_pjm_hs::nhflow_pjm_hs(lexer* p, fdm_nhf *d, patchBC_interface *ppBC)
{
    pBC = ppBC;
    
	pd = new density_f(p);

    
    gcval_press=540;  

	gcval_u=7;
	gcval_v=8;
	gcval_w=9;
    
}

nhflow_pjm_hs::~nhflow_pjm_hs()
{
}

void nhflow_pjm_hs::start(lexer*p, fdm_nhf *d, solver* psolv, ghostcell* pgc, ioflow *pflow, double *U, double *V, double *W, double alpha)
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
    if(p->A521==1 && p->A540==1)
    LOOP
    WETDRY
    d->F[IJK] -= PORVALNH*fabs(p->W22)*
                (p->A523*eta(i+1,j) + (1.0-p->A523)*eta_n(i+1,j) - p->A523*eta(i-1,j) - (1.0-p->A523)*eta_n(i-1,j))/(p->DXP[IP]+p->DXP[IM1]);

    if(p->A521==1 && p->A540==2)
    LOOP
    WETDRY
	d->F[IJK] -= PORVALNH*fabs(p->W22)*(d->eta(i+1,j) - d->eta(i,j))/p->DXP[IP];
    
               
    if(p->A521==2 && p->A540==1)
    LOOP
    WETDRY
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
               
}

void nhflow_pjm_hs::vpgrad(lexer*p, fdm_nhf *d, slice &eta, slice &eta_n)
{
    if(p->A521==1 && p->A540==1)
    LOOP
    WETDRY
	d->G[IJK] -= PORVALNH*fabs(p->W22)*
                 (p->A523*eta(i,j+1) + (1.0-p->A523)*eta_n(i,j+1) - p->A523*eta(i,j-1) - (1.0-p->A523)*eta_n(i,j-1))/(p->DYP[JP]+p->DYP[JM1]);
    
    if(p->A521==1 && p->A540==2)
    LOOP
    WETDRY
	d->G[IJK] -= PORVALNH*fabs(p->W22)*(d->eta(i,j+1) - d->eta(i,j))/p->DYP[JP];
    
    
    if(p->A521==2 && p->A540==1)
    LOOP
    WETDRY
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
}

void nhflow_pjm_hs::wpgrad(lexer*p, fdm_nhf *d, slice &eta, slice &eta_n)
{
}


double nhflow_pjm_hs::limiter(double v1, double v2)
{
    denom = fabs(v1) + fabs(v2);
    
    denom = fabs(denom)>1.0e-10?denom:1.0e10;
    
    val =  (v1*fabs(v2) + fabs(v1)*v2)/denom;

    return val;	
}
