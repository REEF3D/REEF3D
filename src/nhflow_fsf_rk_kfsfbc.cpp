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

#include"nhflow_fsf_rk.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

void nhflow_fsf_rk::kinematic_fsf(lexer *p, fdm_nhf *d, double *U, double *V, double *W, slice &eta1, slice &eta2, double alpha)
{
    double wval,w_n,udetax;
    double Pval,Qval;
    double detax;
    double uvel1,uvel2;
    double zloc1,zloc2;
    
    GC4LOOP
    if(p->gcb4[n][3]==6 && p->gcb4[n][4]==3)
    {
    i=p->gcb4[n][0];
    j=p->gcb4[n][1];
    k=p->gcb4[n][2];
    
    wval = 0.0;
    
        if(p->A515==1)
        {
        Pval = U[IJK];
        Qval = V[IJK];
        
        wval = d->detadt(i,j)
        
             + MAX(0.0,Pval)*((eta1(i,j)-eta1(i-1,j))/(p->DXP[IP]))
             + MIN(0.0,Pval)*((eta1(i+1,j)-eta1(i,j))/(p->DXP[IP1]))
        
             + MAX(0.0,Qval)*((eta1(i,j)-eta1(i,j-1))/(p->DYP[JP]))
             + MIN(0.0,Qval)*((eta1(i,j+1)-eta1(i,j))/(p->DYP[JP1]));
        }
             
        if(p->A515==2)
        {
        wval = d->detadt(i,j)
        
             + U[IJK]*(eta1(i+1,j)-eta1(i-1,j))/(p->DXP[IP]+p->DXP[IP1])

             + V[IJK]*(eta1(i,j+1)-eta1(i,j-1))/(p->DYP[JP]+p->DYP[JP1]);
        }
        
        if(p->A515==3)
        {
        dfdx_plus = (eta1(i+1,j)-eta1(i,j))/p->DXP[IP];
        dfdx_min  = (eta1(i,j)-eta1(i-1,j))/p->DXP[IM1];
    
        detadx = limiter(dfdx_plus,dfdx_min);
        
        dfdy_plus = (eta1(i,j-1)-eta1(i,j))/p->DYP[JP];
        dfdy_min  = (eta1(i,j)-eta1(i,j-1))/p->DYP[JM1];
    
        detady = limiter(dfdy_plus,dfdy_min);
        
        
        wval = d->detadt(i,j)
        
             + U[IJK]*detadx

             + V[IJK]*detady;
        }
        
        //W[IJK] = wval; 
        W[IJKp1] = wval;
        W[IJKp2] = wval;
        W[IJKp3] = wval;
    }
    
    
    
    // Kinematic Bed BC
    GC4LOOP
    if(p->gcb4[n][3]==5 && p->gcb4[n][4]==21)
    {
    i=p->gcb4[n][0];
    j=p->gcb4[n][1];
    k=p->gcb4[n][2];
    
    Pval = U[IJK];
    Qval = V[IJK];
    
    
    if(p->A518==1)
    {
    wval = - MAX(0.0,Pval)*((d->depth(i,j)-d->depth(i-1,j))/(p->DXP[IP]))
           - MIN(0.0,Pval)*((d->depth(i+1,j)-d->depth(i,j))/(p->DXP[IP1]))
           
           - MAX(0.0,Qval)*((d->depth(i,j)-d->depth(i,j-1))/(p->DYP[JP]))
           - MIN(0.0,Qval)*((d->depth(i,j+1)-d->depth(i,j))/(p->DYP[JP1]));
    }
    
    
    if(p->A518==2)
    {
    wval = - Pval*((d->depth(i+1,j)-d->depth(i-1,j))/(p->DXP[IP]+p->DXP[IP1]))
           - Qval*((d->depth(i,j+1)-d->depth(i,j-1))/(p->DYP[JP]+p->DYP[JP1]));
    }
    
        
        //W[IJK] = wval;
        W[IJKm1] = wval;
        W[IJKm2] = wval;
        W[IJKm3] = wval;
        
        w_n = d->wbed(i,j);
        d->wbed(i,j) = wval;
        
        d->dwdt(i,j) = (wval - w_n)/(alpha*p->dt);
    }
}

double nhflow_fsf_rk::limiter(double v1, double v2)
{
    denom = fabs(v1) + fabs(v2);
    
    denom = fabs(denom)>1.0e-10?denom:1.0e10;
    
    val =  (v1*fabs(v2) + fabs(v1)*v2)/denom;

    return val;	
}