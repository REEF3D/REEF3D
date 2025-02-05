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

void grid::fillgcb3(lexer *p)
{
    int q,n;
	
	p->Iarray(p->fgc,imax*jmax*kmax,6);

//  ------------
	
	if(p->gcb3_count!=p->gcb4_count)
	{
	p->Iresize(p->gcb3,p->gcb3_count, p->gcb4_count, 6, 6); 
	p->Dresize(p->gcd3,p->gcb3_count, p->gcb4_count);
	
	p->gcb3_count=p->gcb4_count;
	}
	
	QGCB4
	{
	for(n=0;n<5;++n)
	p->gcb3[q][n]=p->gcb4[q][n];

	if(p->gcb3[q][3]==5 || p->gcb3[q][3]==6)
	p->gcd3[q]=p->gcd4[q];

	if(p->gcb3[q][3]!=5 && p->gcb3[q][3]!=6)
	p->gcd3[q]=p->gcd4[q];
	}

    QGC3LOOP
	{
	    i=p->gcb3[q][0];
		j=p->gcb3[q][1];
		k=p->gcb3[q][2];
        
        if(p->gcb3[q][3]==5 || p->gcb3[q][3]==6)
        p->gcd3[q] += 0.5*p->DZP[KP];

		p->fgc[IJK][p->gcb3[q][3]-1]=1;
	}
	
	
	
	QGC3LOOP
	{
	    i=p->gcb3[q][0];
		j=p->gcb3[q][1];
		k=p->gcb3[q][2];

            if(p->gcb3[q][3]==6 && (p->periodic3!=1 || k+p->origin_k<p->gknoz-1))
            p->gcb3[q][2]-=1;
	}
    
    QGC3LOOP
	{
	    i=p->gcb3[q][0];
		j=p->gcb3[q][1];
		k=p->gcb3[q][2];

            if(p->gcb3[q][3]!=6 && p->fgc[IJK][5]==1 && (p->periodic3!=1 || k+p->origin_k<p->gknoz-1))
            p->gcb3[q][3]=-fabs(p->gcb3[q][3]);
	}
}