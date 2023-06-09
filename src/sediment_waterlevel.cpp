/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"sediment_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment_fdm.h"
#include"topo_relax.h"
#include"ioflow.h"
#include"vrans_v.h"
#include"vrans_f.h"

void sediment_f::waterlevel(lexer *p, fdm *a, ghostcell *pgc)
{
    double zval=0.0;

    SLICELOOP4
    {
        zval = s->bedzh(i,j);

        KLOOP
        PCHECK
        if(a->phi(i,j,k)>=0.0 && a->phi(i,j,k+1)<0.0)
        {
        zval=MAX(zval,-(a->phi(i-1,j,k)*p->DXM)/(a->phi(i,j,k+1)-a->phi(i,j,k)) + p->pos_z());
        }

    s->waterlevel(i,j) = zval-s->bedzh(i,j);
    }
}