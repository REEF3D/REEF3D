/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"nhflow_forcing.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

void nhflow_forcing::solid_forcing(lexer *p, fdm_nhf *d, ghostcell *pgc, 
                             double alpha, double *U, double *V, double *W, slice &WL)
{

// update Heaviside
    pgc->start5V(p,d->FHB,1);
    
    uf=vf=wf=0.0;
    
    LOOP
    {
        H = Hsolidface(p,d,0,0,0);
        
        
        efc = 0.0;
        
        if(d->SOLID[IJK]<0.0)
        {
            efc = 0.0;
            
            if(d->SOLID[Im1JK]>0.0)   
            efc+=1.0;
            
            if(d->SOLID[Ip1JK]>0.0)    
            efc+=1.0;

            if(d->SOLID[IJm1K]>0.0 && p->j_dir==1) 
            efc+=1.0;
            
            if(d->SOLID[IJp1K]>0.0 && p->j_dir==1)    
            efc+=1.0;
            
            if(d->SOLID[IJKm1]>0.0)   
            efc+=1.0;
            
            if(d->SOLID[IJKp1]>0.0)    
            efc+=1.0;
        }
        
        //if(efc>0.1)
        //{
        FX[IJK] += H*(uf - U[IJK])/(alpha*p->dt);
        FY[IJK] += H*(vf - V[IJK])/(alpha*p->dt);
        FZ[IJK] += H*(wf - W[IJK])/(alpha*p->dt);
        
        d->FHB[IJK] = min(d->FHB[IJK] + H, 1.0); 
        //}
    }
    
    pgc->start5V(p,d->FHB,50);
    
    
    k=p->knoz-1;
    /*
    SLICELOOP4
    {
    H = Hsolidface(p,d,0,0,0);
    
    ef = d->depth(i,j);
    
    
    if(d->SOLID[IJK]<0.0)
    {
        ef = 0.0;
        efc = 0.0;
        
        if(d->SOLID[Im1JK]>0.0)   
        {
        ef += WL(i-1,j);
        efc+=1.0;
        }
        
        if(d->SOLID[Ip1JK]>0.0)    
        {
        ef += WL(i+1,j);
        efc+=1.0;
        }

        if(d->SOLID[IJm1K]>0.0 && p->j_dir==1) 
        {
        ef += WL(i,j-1);
        efc+=1.0;
        }
        
        if(d->SOLID[IJp1K]>0.0 && p->j_dir==1)    
        {
        ef += WL(i,j+1);
        efc+=1.0;
        }
        
    if(efc>0.1)
    ef = ef/efc;
    
    if(efc<0.1)
    ef = d->depth(i,j);
    
    }
    
    //if(efc>0.1 && d->SOLID[IJK]<0.0)
    fe(i,j) += H*(ef - WL(i,j))/(alpha*p->dt);
    }*/
}