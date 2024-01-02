/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"nhflow_fsf_f.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

void nhflow_fsf_f::kinematic_fsf(lexer *p, fdm_nhf *d, double *U, double *V, double *W, slice &eta)
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
    
    Pval = d->Ub[IJK];
    Qval = d->Vb[IJK];
    
        if(p->A515==1)
        {
        wval = d->detadt(i,j)
        
             + MAX(0.0,Pval)*((eta(i,j)-eta(i-1,j))/(p->DXP[IP]))
             + MIN(0.0,Pval)*((eta(i+1,j)-eta(i,j))/(p->DXP[IP1]))
        
             + MAX(0.0,Qval)*((eta(i,j)-eta(i,j-1))/(p->DYP[JP]))
             + MIN(0.0,Qval)*((eta(i,j+1)-eta(i,j))/(p->DYP[JP1]));
        }
             
        if(p->A515==2)
        {
        wval = d->detadt(i,j)
        
             + Pval*(eta(i+1,j)-eta(i-1,j))/(p->DXP[IP]+p->DXP[IP1])

             + Qval*(eta(i,j+1)-eta(i,j-1))/(p->DYP[JP]+p->DYP[JP1]);
        }
        
        if(p->A515>=3)
        {
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
        
        d->Wt[IJK] = wval;
        d->Wt[IJKp1] = wval;
        d->Wt[IJKp2] = wval;
        
        if(p->A515==4)
        {
        d->W[IJK] = wval;
        d->W[IJKp1] = wval;
        d->W[IJKp2] = wval;
        }
    }
}   

double nhflow_fsf_f::limiter(double v1, double v2)
{
    denom = fabs(v1) + fabs(v2);
    
    denom = fabs(denom)>1.0e-10?denom:1.0e10;
    
    val =  (v1*fabs(v2) + fabs(v1)*v2)/denom;

    return val;	
}