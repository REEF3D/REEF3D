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

#include"nhflow_sigma.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

#define WLVL (fabs(d->WL(i,j))>(p->A544)?d->WL(i,j):1.0e20)

void nhflow_sigma::sigma_ini(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &eta)
{	

    d->wd_criterion=p->A544;
    

    FLOOP
    p->sig[FIJK] =  p->ZN[KP];
    
    // bc
    SLICELOOP4
    {
        k=0;
        if(p->nb5==-2)
        {
            p->sig[FIJKm1] = p->ZN[KM1];
            p->sig[FIJKm2] = p->ZN[KM2];
            p->sig[FIJKm3] = p->ZN[KM3];
        }
        
        k=p->knoz;
        if(p->nb6==-2)
        {
            p->sig[FIJKp1] = p->ZN[KP1];
            p->sig[FIJKp2] = p->ZN[KP2];
            p->sig[FIJKp3] = p->ZN[KP3];
        } 
    }
    
    pgc->start7S(p,p->sig,1);

    
    SLICELOOP4
    {
    d->Bx(i,j) = 0.0;
    d->By(i,j) = 0.0;
    }
    
    pgc->gcsl_start4(p,d->Bx,50);
    pgc->gcsl_start4(p,d->By,50);
    
    SLICELOOP4
    p->sigz[IJ] = 0.0;
    
    SLICELOOP4
    p->sigt[FIJK] = 0.0;

}

