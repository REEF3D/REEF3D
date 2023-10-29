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

#include"nhflow_fsf_f.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

void nhflow_fsf_f::kinematic_bed(lexer *p, fdm_nhf *d, double *U, double *V, double *W)
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
    
    Pval = U[IJK];
    Qval = V[IJK];
    
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
    
    if(p->A516>=3)
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
    
    if(p->wet[IJ]==0)
    wval=0.0;

        d->Wb[IJKm1] = wval;
        d->Wb[IJKm2] = wval;
        d->Wb[IJKm3] = wval;

        d->Wt[IJKm1] = wval;
        d->Wt[IJKm2] = wval;
        d->Wt[IJKm3] = wval;
        
        if(p->A516==4)
        {
        d->W[IJKm1] = wval;
        d->W[IJKm2] = wval;
        d->W[IJKm3] = wval;
        }
        
    }
}

