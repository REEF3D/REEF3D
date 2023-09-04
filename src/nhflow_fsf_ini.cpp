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
#include"ioflow.h"

void nhflow_fsf_f::ini(lexer* p, fdm_nhf* d, ghostcell* pgc, ioflow* pflow, double *U, double *V, double *W)
{   
    if(d->WL(i,j)<=p->A544)
    {
        p->wet[IJ]=0;
        d->eta(i,j) =  p->A544 - d->depth(i,j) - eps;
        d->WL(i,j) = d->eta(i,j) + d->depth(i,j);
    }
    
    pgc->gcsl_start4Vint(p,p->wet,50);
    
    wetdry(p,d,pgc,U,V,W,d->eta);
    
    SLICELOOP4
    d->detadt(i,j) = 0.0;
    
    pgc->gcsl_start4(p,d->detadt,1);
    pgc->start1V(p,d->Fx,10);
    pgc->start2V(p,d->Fy,10);
    
    LOOP
    d->detadt(i,j) += -p->DZN[KP]*((d->Fx[IJK] - d->Fx[Im1JK])/p->DXN[IP]  + (d->Fy[IJK] - d->Fy[IJm1K])/p->DYN[JP]*p->y_dir);
    
    pgc->gcsl_start4(p,d->detadt,1);
    
    pgc->start1V(p,d->Fx,10);
    
    //LOOP    
    //d->test[IJK] = d->Fx[IJK] - d->Fx[Im1JK];
    
    pgc->start4V(p,d->test,1);
     
}