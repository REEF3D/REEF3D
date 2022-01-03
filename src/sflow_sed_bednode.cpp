/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2022 Hans Bihs

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

#include"sflow_sediment_f.h"
#include"lexer.h"
#include"fdm2D.h" 
#include"ghostcell.h"

void sflow_sediment_f::bednode(lexer *p, fdm2D *b, ghostcell *pgc)
{
    pgc->gcsl_start4(p,b->bed,50);

    // bednode upate
    TPSLICELOOP
    {
    pip=4;
    
    b->bednode(i,j) = 0.25*(b->bed(i,j) + b->bed(i+1,j) + b->bed(i,j+1) + b->bed(i,j));

    if(p->flagslice4[Im1Jm1]<0 && p->flagslice4[Im1J]>0 && p->flagslice4[IJm1]>0 && p->flagslice4[IJ]>0)
    b->bednode(i,j) = (1.0/3.0)*(b->bed(i+1,j) + b->bed(i,j+1) + b->bed(i,j));
    
    if(p->flagslice4[Im1Jm1]>0 && p->flagslice4[Im1J]<0 && p->flagslice4[IJm1]>0 && p->flagslice4[IJ]>0)
    b->bednode(i,j) = (1.0/3.0)*(b->bed(i+1,j+1) + b->bed(i,j+1) + b->bed(i,j));
    
    if(p->flagslice4[Im1Jm1]>0 && p->flagslice4[Im1J]>0 && p->flagslice4[IJm1]<0 && p->flagslice4[IJ]>0)
    b->bednode(i,j) = (1.0/3.0)*(b->bed(i+1,j+1) + b->bed(i+1,j) + b->bed(i,j));
    
    if(p->flagslice4[Im1Jm1]>0 && p->flagslice4[Im1J]>0 && p->flagslice4[IJm1]>0 && p->flagslice4[IJ]<0)
    b->bednode(i,j) = (1.0/3.0)*(b->bed(i+1,j+1) + b->bed(i+1,j) + b->bed(i,j+1));
    
    pip=0;
    } 
    
    
}
