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
along with this program; if not, see <http://www.gnu.org/liceonephases/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"ptf_bed_update.h"
#include"fnpf_fsfbc.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"fnpf_cds2.h"
#include"fnpf_cds4.h"

ptf_bed_update::ptf_bed_update(lexer *p, fdm *a, ghostcell *pgc) 
{    
    if(p->A313==2)
    pconvec = new fnpf_cds2(p);
    
    if(p->A313==3)
    pconvec = new fnpf_cds4(p);
}

ptf_bed_update::~ptf_bed_update()
{
}

void ptf_bed_update::bedbc(lexer *p, fdm *a, ghostcell *pgc, field &Fi)
{
    double Fval;
    
    GC4LOOP
    if(p->gcb4[n][3]==5 && (p->gcb4[n][4]==5 || p->gcb4[n][4]==21))
    {
        i = p->gcb4[n][0];
        j = p->gcb4[n][1];
        k = p->gcb4[n][2];
        
        Fval = ((a->depth(i+1,j)-a->depth(i-1,j))/(p->DXP[IP]+p->DXP[IM1]))*((Fi(i+1,j,k)-Fi(i-1,j,k))/(p->DXP[IP]+p->DXP[IM1]))
             + ((a->depth(i,j+1)-a->depth(i,j-1))/(p->DYP[JP]+p->DYP[JM1]))*((Fi(i,j+1,k)-Fi(i,j-1,k))/(p->DYP[JP]+p->DYP[JM1]));
    

        Fi(i,j,k-1) = Fi(i,j,k) + 0.0*p->DXM*Fval;
        Fi(i,j,k-2) = Fi(i,j,k) + 0.0*p->DXM*Fval;
        Fi(i,j,k-3) = Fi(i,j,k) + 0.0*p->DXM*Fval;
    }
}


void ptf_bed_update::waterdepth(lexer *p, fdm *a, ghostcell *pgc)
{
    SLICELOOP4
	a->depth(i,j) = p->wd - a->bed(i,j);
}




