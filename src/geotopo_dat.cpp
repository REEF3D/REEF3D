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

#include"lexer.h"
#include"geotopo.h"
#include"fdm.h"
#include"ghostcell.h"

void geotopo::dat(lexer* p, fdm* a, ghostcell* pgc)
{
    if(p->toporead>0)
    BASELOOP
    a->topo(i,j,k) = p->flag_topo[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin];
    
    if(p->S57>-1.0e20)
    ALOOP
    a->topo(i,j,k)=-p->S57+p->ZP[KP];
    
    pgc->start4a(p,a->topo,150);
    
    p->del_Darray(p->flag_topo,p->imax*p->jmax*p->kmax);
    
}






