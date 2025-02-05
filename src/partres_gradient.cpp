/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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
Authors: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#include"partres.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment_fdm.h"

void partres::stress_gradient(lexer *p, fdm *a, ghostcell *pgc, sediment_fdm *s)
{
    ALOOP
    {
    dTx(i,j,k) = ((Tau(i+1,j,k) - Tau(i-1,j,k))/(p->DXP[IM1]+p->DXP[IP]))/((Tsval>1.0e-6?Tsval:1.0e10));
    dTy(i,j,k) = ((Tau(i,j+1,k) - Tau(i,j-1,k))/(p->DYP[JM1]+p->DYP[JP]))/((Tsval>1.0e-6?Tsval:1.0e10));
    dTz(i,j,k) = ((Tau(i,j,k+1) - Tau(i,j,k-1))/(p->DZP[KM1]+p->DZP[KP]))/((Tsval>1.0e-6?Tsval:1.0e10));
    }
    
    pgc->start4a(p,dTx,1);
    pgc->start4a(p,dTy,1);
    pgc->start4a(p,dTz,1);
}

void partres::pressure_gradient(lexer *p, fdm *a, ghostcell *pgc, sediment_fdm *s)
{
    ALOOP
    {
    dPx(i,j,k) = (a->press(i+1,j,k) - a->press(i-1,j,k))/(p->DXP[IM1]+p->DXP[IP]);
    dPy(i,j,k) = (a->press(i,j+1,k) - a->press(i,j-1,k))/(p->DYP[JM1]+p->DYP[JP]);
    dPz(i,j,k) = (a->press(i,j,k+1) - a->press(i,j,k-1))/(p->DZP[KM1]+p->DZP[KP]);
    }
    
    pgc->start4a(p,dPx,1);
    pgc->start4a(p,dPy,1);
    pgc->start4a(p,dPz,1);
}

