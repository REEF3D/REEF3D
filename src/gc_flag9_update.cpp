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

#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"

void ghostcell::flag9_update(lexer *p, fdm *a)
{
    
    for(i=0;i<p->imax*p->jmax*p->kmax; ++i)
    p->flag9[i]=1;
    
    BASELOOP
    if(a->solid(i,j,k)<0.0 || a->topo(i,j,k)<0.0)
    p->flag9[IJK]=-1;
    
    
    
    LOOP
    {
    if(a->phi(i,j,k)>0 && (a->solid(i,j,k-1)<0.0 || a->topo(i,j,k-1)<0.0))
    {
    p->flag9[IJKm1]=1;
    p->flag9[IJKm2]=1;
    p->flag9[IJKm3]=1;
    }
    
    if(a->phi(i,j,k)>0 && a->phi(i,j,k-1)>0 && (a->solid(i,j,k-2)<0.0 || a->topo(i,j,k-2)<0.0))
    {
    p->flag9[IJKm2]=1;
    p->flag9[IJKm3]=1;
    }
        
        
        
    }
    
    sizeM_update(p,a);
    
}