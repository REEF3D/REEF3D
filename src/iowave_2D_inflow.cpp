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

#include"iowave.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"patchBC_interface.h"

void iowave::inflow2D(lexer *p, fdm2D* b, ghostcell* pgc, slice &P, slice &Q, slice &bed, slice &eta)
{
    if(p->B98==0)
    inflow2D_plain(p,b,pgc,P,Q,eta);
    
	if(p->B98==3 || p->B98==4)
	wavegen2D(p,b,pgc,P,Q,bed,eta);
	
	if(p->B99==3 || p->B99==4)
	active_beach2D(p,b,pgc,P,Q,bed,eta);
    
    pBC->patchBC_ioflow2D(p,pgc,P,Q,bed,eta);
}

void iowave::rkinflow2D(lexer *p, fdm2D* b, ghostcell* pgc, slice &P, slice &Q, slice &U, slice &V)
{
    for(n=0;n<p->gcslin_count;n++)
    {
    i=p->gcslin[n][0];
    j=p->gcslin[n][1];


        P(i-1,j)=U(i-1,j);
        P(i-2,j)=U(i-2,j);
        P(i-3,j)=U(i-3,j);

        Q(i-1,j)=V(i-1,j);
        Q(i-2,j)=V(i-2,j);
        Q(i-3,j)=V(i-3,j);
    }
    
    pBC->patchBC_rkioflow2D(p,pgc,P,Q,U,V);
}

void iowave::inflow2D_plain(lexer *p, fdm2D* b, ghostcell* pgc, slice &P, slice &Q, slice &eta)
{
    for(n=0;n<p->gcslin_count;n++)
    {
    i=p->gcslin[n][0];
    j=p->gcslin[n][1];

        P(i-1,j)=p->Ui;
        P(i-2,j)=p->Ui;
        P(i-3,j)=p->Ui;

		Q(i-1,j)=0.0;
        Q(i-2,j)=0.0;
        Q(i-3,j)=0.0;
        
        eta(i-1,j)=0.0;
        eta(i-2,j)=0.0;
        eta(i-3,j)=0.0;

        b->hx(i-1,j)=b->eta(i-1,j) + b->depth(i-1,j);
        b->hx(i-2,j)=b->eta(i-2,j) + b->depth(i-2,j);
        b->hx(i-3,j)=b->eta(i-3,j) + b->depth(i-3,j);
        
        b->hy(i-1,j)=b->eta(i-1,j) + b->depth(i-1,j);
        b->hy(i-2,j)=b->eta(i-2,j) + b->depth(i-2,j);
        b->hy(i-3,j)=b->eta(i-3,j) + b->depth(i-3,j);
    }
}
