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

#include"nhflow_force.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include <math.h>

void nhflow_force::allocate(lexer* p, fdm_nhf *d, ghostcell *pgc)
{
    p->Iarray(tri,numtri,4);
    p->Darray(pt,numvert,3);
    p->Darray(ls,numvert);
    p->Iarray(facet,numtri,4);
    p->Iarray(confac,numtri);
    p->Iarray(numfac,numtri);
	p->Iarray(numpt,numtri);
    p->Darray(ccpt,numtri*4,3);
    
    
    // ini
    int n,q;
    
    for(n=0; n<numtri; ++n)
    {
    for(q=0;q<4;++q)
    tri[n][q]=0;
    
    for(q=0;q<4;++q)
    facet[n][q]=0;
    
    
    confac[n]=0;
    numfac[n]=0;
    numpt[n]=0;
    }
    
    for(n=0; n<numvert; ++n)
    {
    for(q=0;q<3;++q)
    pt[n][q]=0.0;
    
    ls[n]=0.0;
    }
    
    for(n=0; n<numtri*4; ++n)
    {
    for(q=0;q<3;++q)
    ccpt[n][q]=0;
    }
}

void nhflow_force::deallocate(lexer* p, fdm_nhf *d, ghostcell *pgc)
{
    p->del_Iarray(tri,numtri,4);
    p->del_Darray(pt,numvert,3);
    p->del_Darray(ls,numvert);
    p->del_Iarray(facet,numtri,4);
    p->del_Iarray(confac,numtri);
    p->del_Iarray(numfac,numtri);
	p->del_Iarray(numpt,numtri);
    p->del_Darray(ccpt,numtri*4,3);
}
