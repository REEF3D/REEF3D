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
#include"ioflow.h"
#include"patchBC_interface.h"

void nhflow_fsf_f::ucorr(lexer* p, fdm_nhf* d, double *UH, slice &WL, double alpha)
{
    /*if(p->A521==1)
    LOOP
    {
        dfdx_plus = (d->detadt(i+1,j)-d->detadt(i,j))/p->DXP[IP];
        dfdx_min  = (d->detadt(i,j)-d->detadt(i-1,j))/p->DXP[IM1];
    
        detadx = limiter(dfdx_plus,dfdx_min);
        
    UH[IJK] += 0.25*alpha*alpha*p->dt*p->dt*WL(i,j)*fabs(p->W22)*detadx;
    }*/
}

void nhflow_fsf_f::vcorr(lexer* p, fdm_nhf* d, double *VH, slice &WL, double alpha)
{
    /*if(p->A521==1 && p->j_dir==1)
    LOOP
    {
        dfdy_plus = (d->detadt(i,j+1)-d->detadt(i,j))/p->DYP[JP];
        dfdy_min  = (d->detadt(i,j)-d->detadt(i,j-1))/p->DYP[JM1];
    
        detady = limiter(dfdy_plus,dfdy_min);
        
    VH[IJK] += 0.25*alpha*alpha*p->dt*p->dt*WL(i,j)*fabs(p->W22)*detady;
    }*/
}