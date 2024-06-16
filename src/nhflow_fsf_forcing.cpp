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

void nhflow_fsf_f::forcing(lexer* p, fdm_nhf* d, ghostcell* pgc, double *UH, double *VH, double *WH, slice &WL)
{/*
    if(p->A561>0 || p->A564>0)
    {
     
    k=p->knoz-1;
    
    SLICELOOP4
    if(d->SOLID[IJK]>0)
    {
    
        if(d->SOLID[Im1JK]<0)
        {
        WL(i-1,j) = WL(i,j);
        d->eta(i-1,j) = d->eta(i,j);
        }
        
        if(d->SOLID[Ip1JK]<0)
        {
        WL(i+1,j) = WL(i,j);
        d->eta(i+1,j) = d->eta(i,j);
        }
        
        
    }
    
    pgc->gcsl_start4(p,WL,gcval_eta);
    pgc->gcsl_start4(p,d->eta,gcval_eta);

    }*/
    
}
