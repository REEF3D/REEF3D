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
Authors: Hans Bihs
--------------------------------------------------------------------*/

#include"CPM.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment_fdm.h"

void CPM::press_lithostatic(lexer *p, fdm *a, ghostcell *pgc, sediment_fdm *s)
{
    ALOOP
    press(i,j,k) = 0.0;
    
    pgc->start4a(p,press,1);
    
    ILOOP
    JLOOP
    {    
        KREVLOOP
        {
        //if(a->topo(i,j,k-2)<0.0)
        press(i,j,k) = fabs(p->W22) * (p->S22 - a->ro(i,j,k)) * (-a->topo(i,j,k));
            
        }
        

    }

    
    pgc->start4a(p,press,1);
}
