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

#include"solid.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"reinitopo.h"
#include"ioflow.h"
#include"reinitopo_RK3.h"

solid::solid(lexer* p, fdm *a, ghostcell* pgc)
{
}

solid::~solid()
{
}

void solid::start(lexer* p, fdm* a, ghostcell* pgc, ioflow *pflow, convection* pconvec, reinitopo* preso)
{

	solid_topo(p,a,pgc);
    
    preso->start(p,a,pgc,a->solid);
    
    if(p->G3==0)
    {
    pgc->solid_update(p,a);

    pflow->gcio_update(p,a,pgc);
    }
}

void solid::solid_topo(lexer* p, fdm* a, ghostcell* pgc)
{
    BASELOOP
    {
    a->solid(i,j,k) = p->flag_solid[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin];
    }
    
    p->del_Darray(p->flag_solid,p->imax*p->jmax*p->kmax);
}

