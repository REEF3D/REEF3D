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

#define WLVL (fabs(d->WL(i,j))>1.0e-20?d->WL(i,j):1.0-20)

void nhflow_sigma::sigma_ini(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &eta)
{	
    SLICELOOP4
    p->wet[IJ]=1;
    
    pgc->gcsl_start4Vint(p,p->wet,50);
    
    SLICELOOP4
    p->wet_n[IJ]=1;
    
    pgc->gcsl_start4Vint(p,p->wet_n,50);
    
    
    d->wd_criterion=0.00005;
    
    p->Darray(p->sig, p->imax*p->jmax*(p->kmax+2));
    p->Darray(p->sigx,p->imax*p->jmax*(p->kmax+2));
    p->Darray(p->sigy,p->imax*p->jmax*(p->kmax+2));
    p->Darray(p->sigz,p->imax*p->jmax);
    p->Darray(p->sigt,p->imax*p->jmax*(p->kmax+2));

    p->Darray(p->sigxx,p->imax*p->jmax*(p->kmax+2));
    

    FLOOP
    p->sig[FIJK] =  p->ZN[KP];
    
    // bc
    SLICELOOP4
    {
        k=0;
        if(p->nb5==-2)
        {
            p->sig[FIJKm1] = p->ZN[FIJK];
            p->sig[FIJKm2] = p->ZN[FIJK];
            p->sig[FIJKm3] = p->ZN[FIJK];
        }
        
        k=p->knoz;
        if(p->nb6==-2)
        {
            p->sig[FIJKp1] = p->ZN[FIJK];
            p->sig[FIJKp2] = p->ZN[FIJK];
            p->sig[FIJKp3] = p->ZN[FIJK];
        } 
    }
    
    pgc->start7S(p,p->sig,1);

    
    SLICELOOP4
	d->bed(i,j) = p->bed[IJ];
    
    SLICELOOP4
    d->WL(i,j) = MAX(p->A544, d->eta(i,j) + p->wd - d->bed(i,j));
    
    //SLICELOOP4
	//d->depth(i,j) = MAX(p->wd - d->bed(i,j), p->A544);
    
    SLICELOOP4
	d->depth(i,j) = p->wd - d->bed(i,j);
    
    SLICELOOP4
    {
    d->Bx(i,j) = 0.0;
    d->By(i,j) = 0.0;
    }
    
    pgc->gcsl_start4(p,d->depth,50);
    pgc->gcsl_start4(p,d->Bx,50);
    pgc->gcsl_start4(p,d->By,50);
    
    SLICELOOP4
    p->sigz[IJ] = 1.0/WLVL;
    
    SLICELOOP4
    p->sigt[FIJK] = 0.0;

}


