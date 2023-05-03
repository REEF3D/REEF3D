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
along with this program; if not, see <http://www.gnu.org/liceonephases/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"fnpf_bed_update.h"
#include"fnpf_fsf.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"fnpf_cds2.h"
#include"fnpf_cds4.h"

fnpf_bed_update::fnpf_bed_update(lexer *p) 
{    
    if(p->A313==2)
    pconvec = new fnpf_cds2(p);
    
    if(p->A313==3)
    pconvec = new fnpf_cds4(p);
}

fnpf_bed_update::~fnpf_bed_update()
{
}

void fnpf_bed_update::bedbc_sig(lexer *p, fdm_fnpf *c, ghostcell *pgc, double *Fi, fnpf_fsf *pf)
{
    // Calculate gradients
    
    double bcval,denom;
    k=0;
    SLICELOOP4
    {
    bcval =   (c->Bx(i,j)*(c->Fi[FIp1JK]-c->Fi[FIm1JK])/(p->DXP[IP] + p->DXP[IM1])
    
            +  c->By(i,j)*(c->Fi[FIJp1K]-c->Fi[FIJm1K])/(p->DYP[JP] + p->DYP[JM1]));
            
    denom =  p->sigz[IJ] + c->Bx(i,j)*p->sigx[FIJK] + c->By(i,j)*p->sigy[FIJK];
    
    Fi[FIJKm1] = (bcval/denom)*(1.0*p->DZN[KP]) + Fi[FIJK];
    Fi[FIJKm2] = (bcval/denom)*(2.0*p->DZN[KP]) + Fi[FIJK];
    Fi[FIJKm3] = (bcval/denom)*(3.0*p->DZN[KP]) + Fi[FIJK];
    }
}

void fnpf_bed_update::waterdepth(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    SLICELOOP4
	c->depth(i,j) = p->wd - c->bed(i,j);
}




