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

#include"nhflow_sigma.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

void nhflow_sigma::bed_acceleration(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &WL, double *U, double *V, double *W)
{ 
    
    double du,dv;
    k=0;
    
    
    SLICELOOP4
    {
    du = (U[IJK] - d->un(i,j))/p->dt;
    dv = (V[IJK] - d->vn(i,j))/p->dt;
    
    d->dudt(i,j) = sqrt(du*du + dv*dv);
    }
        
    SLICELOOP4
    {
    d->un(i,j) = U[IJK];
    d->vn(i,j) = V[IJK];
    }
    
    
    
    /*
    double du,dv;
    double Uda,Vda;
    k=0;
    SLICELOOP4
    {
        
        Uda=Vda=0.0;
        KLOOP
        {
        Uda += d->U[IJK]*p->DZN[KP]*d->WL(i,j);
        Vda += d->V[IJK]*p->DZN[KP]*d->WL(i,j);
        }
        
        Uda=Uda/d->WL(i,j);
        Vda=Vda/d->WL(i,j);
            
    du = (Uda - d->un(i,j))/p->dt;
    dv = (Vda - d->vn(i,j))/p->dt;
    
    d->dudt(i,j) = sqrt(du*du + dv*dv);
   
    d->un(i,j) = Uda;
    d->vn(i,j) = Vda;
    }*/
    
    
}