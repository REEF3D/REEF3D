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

void nhflow_fsf_rk::kinematic_fsf(lexer *p, fdm_nhf *d, double *U, double *V, double *W, slice &eta)
{
    double wval,w_n;
    double Pval,Qval;
    
    GC4LOOP
    if(p->gcb4[n][3]==6 && p->gcb4[n][4]==3)
    {
    i=p->gcb4[n][0];
    j=p->gcb4[n][1];
    k=p->gcb4[n][2];
    
    wval = 0.0;
    
        if(p->A515==1)
        {
        Pval = d->Ub[IJK];
        Qval = d->Vb[IJK];
        
        wval = d->detadt(i,j)
        
             + MAX(0.0,Pval)*((eta(i,j)-eta(i-1,j))/(p->DXP[IP]))
             + MIN(0.0,Pval)*((eta(i+1,j)-eta(i,j))/(p->DXP[IP1]))
        
             + MAX(0.0,Qval)*((eta(i,j)-eta(i,j-1))/(p->DYP[JP]))
             + MIN(0.0,Qval)*((eta(i,j+1)-eta(i,j))/(p->DYP[JP1]));
        }
             
        if(p->A515==2)
        {
        Pval = d->Ub[IJK];
        Qval = d->Vb[IJK];
        
        wval = d->detadt(i,j)
        
             + Pval*(eta(i+1,j)-eta(i-1,j))/(p->DXP[IP]+p->DXP[IP1])

             + Qval*(eta(i,j+1)-eta(i,j-1))/(p->DYP[JP]+p->DYP[JP1]);
        }
        
        if(p->A515==3)
        {
        Pval = d->Ub[IJK];
        Qval = d->Vb[IJK];
        
        dfdx_plus = (eta(i+1,j)-eta(i,j))/p->DXP[IP];
        dfdx_min  = (eta(i,j)-eta(i-1,j))/p->DXP[IM1];
    
        detadx = limiter(dfdx_plus,dfdx_min);
        
        dfdy_plus = (eta(i,j+1)-eta(i,j))/p->DYP[JP];
        dfdy_min  = (eta(i,j)-eta(i,j-1))/p->DYP[JM1];
    
        detady = limiter(dfdy_plus,dfdy_min);
        
        
        wval = d->detadt(i,j)
        
             + Pval*detadx

             + Qval*detady;
        }

        d->Wb[IJK] = wval;
        d->Wb[IJKp1] = wval;
        d->Wb[IJKp2] = wval;
        d->Wb[IJKp3] = wval;
        
        d->Wt[IJKp1] = wval;
        d->Wt[IJKp2] = wval;
        d->Wt[IJKp3] = wval;
    }
}   
  
void nhflow_fsf_rk::kinematic_bed(lexer *p, fdm_nhf *d, double *U, double *V, double *W)
{
    double wval,w_n;
    double Pval,Qval;
  
    // Kinematic Bed BC
    GC4LOOP
    if(p->gcb4[n][3]==5 && p->gcb4[n][4]==21)
    {
    i=p->gcb4[n][0];
    j=p->gcb4[n][1];
    k=p->gcb4[n][2];
    
    Pval = d->Ut[IJK];
    Qval = d->Vt[IJK];
    
    wval=0.0;    
    
    if(p->A516==1)
    {
    wval = - MAX(0.0,Pval)*((d->depth(i,j)-d->depth(i-1,j))/(p->DXP[IP]))
           - MIN(0.0,Pval)*((d->depth(i+1,j)-d->depth(i,j))/(p->DXP[IP1]))
           
           - MAX(0.0,Qval)*((d->depth(i,j)-d->depth(i,j-1))/(p->DYP[JP]))
           - MIN(0.0,Qval)*((d->depth(i,j+1)-d->depth(i,j))/(p->DYP[JP1]));
    }
    
    
    if(p->A516==2)
    {
    wval = - Pval*((d->depth(i+1,j)-d->depth(i-1,j))/(p->DXP[IP]+p->DXP[IP1]))
           - Qval*((d->depth(i,j+1)-d->depth(i,j-1))/(p->DYP[JP]+p->DYP[JP1]));
    }
    
    if(p->A516==3)
    {
        dfdx_plus = (d->depth(i+1,j)-d->depth(i,j))/p->DXP[IP];
        dfdx_min  = (d->depth(i,j)-d->depth(i-1,j))/p->DXP[IM1];
    
        detadx = limiter(dfdx_plus,dfdx_min);
        
        dfdy_plus = (d->depth(i,j+1)-d->depth(i,j))/p->DYP[JP];
        dfdy_min  = (d->depth(i,j)-d->depth(i,j-1))/p->DYP[JM1];
    
        detady = limiter(dfdy_plus,dfdy_min);
        
        
        wval = - Pval*detadx

               - Qval*detady;
    }
    
        
        d->Wb[IJKm1] = wval;
        d->Wb[IJKm2] = wval;
        d->Wb[IJKm3] = wval;

        d->Wt[IJKm1] = wval;
        d->Wt[IJKm2] = wval;
        d->Wt[IJKm3] = wval;
    }
}

double nhflow_fsf_rk::limiter(double v1, double v2)
{
    denom = fabs(v1) + fabs(v2);
    
    denom = fabs(denom)>1.0e-10?denom:1.0e10;
    
    val =  (v1*fabs(v2) + fabs(v1)*v2)/denom;

    return val;	
}