/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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
Author: Alexander Hanke, Hans Bihs
--------------------------------------------------------------------*/

#include"partres.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment_fdm.h"

void partres::cellSum_update(lexer *p, ghostcell *pgc, sediment_fdm *s, double *PX, double *PY, double *PZ)
{
    ALOOP
    cellSum(i,j,k) = 0.0;
    
    pgc->start4a(p,cellSum,1);

    double Sx,Sy,Sz;

    for(size_t n=0;n<P.index;n++)
    if(P.Flag[n]>=ACTIVE)
    {
        i=p->posf_i(P.XRK1[n]);
        j=p->posf_j(P.YRK1[n]);
        k=p->posf_k(P.ZRK1[n]);

        Sx = (p->XP[IP1] - PX[n])/p->DXP[IP];
        Sy = (p->YP[JP1] - PY[n])/p->DYP[JP];
        Sz = (p->ZP[KP1] - PZ[n])/p->DZP[KP];

        cellSum(i,j,k) += P.ParcelFactor * Sx*Sy*Sz;
        cellSum(i+1,j,k) += P.ParcelFactor * (1.0-Sx)*Sy*Sz;
        cellSum(i+1,j+1,k) += P.ParcelFactor * (1.0-Sx)*(1.0-Sy)*Sz;
        cellSum(i,j+1,k) += P.ParcelFactor * Sx*(1.0-Sy)*Sz;
        cellSum(i,j,k+1) += P.ParcelFactor * Sx*Sy*(1.0-Sz);
        cellSum(i+1,j,k+1) += P.ParcelFactor * (1.0-Sx)*Sy*(1.0-Sz);
        cellSum(i+1,j+1,k+1) += P.ParcelFactor * (1.0-Sx)*(1.0-Sy)*(1.0-Sz);
        cellSum(i,j+1,k+1) += P.ParcelFactor * Sx*(1.0-Sy)*(1.0-Sz);
    }
    

    pgc->start4a(p,cellSum,1);
    
    pgc->start4a_sum(p,cellSum,1);
    
    ALOOP
    Ts(i,j,k) = (1.0/6.0)*PI*pow(P.d50,3.0)*cellSum(i,j,k)/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
    
    pgc->start4a(p,Ts,1);
}

