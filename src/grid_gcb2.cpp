/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"grid.h"
#include"lexer.h"
#include"ghostcell.h"

void grid::fillgcb2(lexer *p)
{
    int q,n;
	
	p->Iarray(p->fgc,imax*jmax*kmax,6);
	
// ----

	if(p->gcb2_count!=p->gcb4_count)
	{
	p->Iresize(p->gcb2,p->gcb2_count, p->gcb4_count, 6, 6); 
	p->Dresize(p->gcd2,p->gcb2_count, p->gcb4_count); 
	
	p->gcb2_count=p->gcb4_count;
	}

	QGCB4
	{
	for(n=0;n<5;++n)
	p->gcb2[q][n]=p->gcb4[q][n];

	if(p->gcb2[q][3]==2 || p->gcb2[q][3]==3)
	p->gcd2[q]=p->gcd4[q];

	if(p->gcb2[q][3]!=2 && p->gcb2[q][3]!=3)
	p->gcd2[q]=p->gcd4[q];
	}

    QGC2LOOP
	{
	    i=p->gcb2[q][0];
		j=p->gcb2[q][1];
		k=p->gcb2[q][2];
        
        if(p->gcb2[q][3]==2 || p->gcb2[q][3]==3)
        p->gcd2[q] += 0.5*p->DYP[JP];

		p->fgc[IJK][p->gcb2[q][3]-1]=1;
	}
	
	
	
	QGC2LOOP
	{
	    i=p->gcb2[q][0];
		j=p->gcb2[q][1];
		k=p->gcb2[q][2];

            if(p->gcb2[q][3]==2 && (p->periodic2!=1 || j+p->origin_j<p->gknoy-1))
            p->gcb2[q][1]-=1;
	}
    
    QGC2LOOP
	{
	    i=p->gcb2[q][0];
		j=p->gcb2[q][1];
		k=p->gcb2[q][2];

            if(p->gcb2[q][3]!=2 && p->fgc[IJK][1]==1 && (p->periodic2!=1 || j+p->origin_j<p->gknoy-1))
            p->gcb2[q][3]=-fabs(p->gcb2[q][3]);
	}
}