/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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
    
        // curvilinear grid: inflow along the normal of the inflow face
        if(p->wet[IJ]==1 && p->geo_curv==1)
        {
        const int qf = (i-1-p->imin)*p->jmax + (j-p->jmin);
        const double nx = p->SX1[qf]/p->LEN1[qf], ny = p->SY1[qf]/p->LEN1[qf];
        const double ui = (p->B62==1) ? Uk*pow(MAX(b->hp(i-1,j),0.0),2.0/3.0) : p->Ui;
        
        P(i-1,j)=P(i-2,j)=P(i-3,j)=ui*nx;
        Q(i-1,j)=Q(i-2,j)=Q(i-3,j)=ui*ny;
        continue;
        }
        
        if(p->wet[IJ]==1)
        {
        const double ui = (p->B62==1) ? Uk*pow(MAX(b->hp(i-1,j),0.0),2.0/3.0) : p->Ui;
        P(i-1,j)=ui;
        P(i-2,j)=ui;
        P(i-3,j)=ui;
        }
        
        if(p->wet[IJ]==0)// || p->gcslin[n][5]==0)
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
