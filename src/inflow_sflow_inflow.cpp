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

#include"ioflow_f.h"
#include"lexer.h"
#include"ghostcell.h"
#include"slice.h"
#include"fdm2D.h"
#include"patchBC_interface.h"

void ioflow_f::inflow2D(lexer *p, fdm2D* b, ghostcell* pgc, slice &P, slice &Q, slice &bed, slice &eta)
{
    for(n=0;n<p->gcslin_count;n++)
    {
    i=p->gcslin[n][0];
    j=p->gcslin[n][1];
    
        if(p->wet[IJ]==1 && p->gcslin[n][5]==1)
        {
        P(i-1,j)=p->Ui;
        P(i-2,j)=p->Ui;
        P(i-3,j)=p->Ui;
        }
        
        if(p->wet[IJ]==0 || p->gcslin[n][5]==0)
        {
        P(i-1,j)=0.0;
        P(i-2,j)=0.0;
        P(i-3,j)=0.0;
        }
		
		Q(i-1,j)=0.0;
        Q(i-2,j)=0.0;
        Q(i-3,j)=0.0;
    }
    
    pBC->patchBC_ioflow2D(p,pgc,P,Q,bed,eta);
}

void ioflow_f::rkinflow2D(lexer *p, fdm2D* b, ghostcell* pgc, slice &P, slice &Q, slice &U, slice &V)
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
