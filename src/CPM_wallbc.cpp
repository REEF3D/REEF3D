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

#include"CPM.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment_fdm.h"

// parcels in the first cell column at the in- and outflow boundaries are fixed (BEDBC),
// they retain the bed at the open boundaries
void CPM::wallbc(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
    for(size_t n=0;n<P.index;n++)
    if(P.Flag[n]==ACTIVE)
    {
        i=p->posc_i(P.X[n]);
        
        if((i<=0 && p->nb1<0 && open_side[0]==1) || (i>=p->knox-1 && p->nb4<0 && open_side[3]==1))
        {
        P.Flag[n]=BEDBC;
        
        P.U[n] = P.V[n] = P.W[n] = 0.0;
        
        P.XRK1[n] = P.X[n];
        P.YRK1[n] = P.Y[n];
        P.ZRK1[n] = P.Z[n];
        
        P.URK1[n] = P.VRK1[n] = P.WRK1[n] = 0.0;
        }
    }
}
