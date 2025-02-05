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
Author: Alexander Hanke, Hans Bihs
--------------------------------------------------------------------*/

#include"partres.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment_fdm.h"

void partres::cellSum_update(lexer *p, ghostcell *pgc, sediment_fdm *s, int mode)
{
    for(n=0;n<P.index;++n)
    if(P.Flag[n]==ACTIVE)
    {
        // step 1
        if(mode==1)
        {
        i=p->posc_i(P.X[n]);
        j=p->posc_j(P.Y[n]);
        k=p->posc_k(P.Z[n]);
        }
        
        if(mode==2)
        {
        i=p->posc_i(P.XRK1[n]);
        j=p->posc_j(P.YRK1[n]);
        k=p->posc_k(P.ZRK1[n]);
        }
        
        cellSum(i,j,k) -= P.ParcelFactor;
        bedch(i,j) -= P.ParcelFactor;
        
        
        // step 2
        if(mode==1)
        {
        i=p->posc_i(P.XRK1[n]);
        j=p->posc_j(P.YRK1[n]);
        k=p->posc_k(P.ZRK1[n]);
        }
        
        if(mode==2)
        {
        i=p->posc_i(P.X[n]);
        j=p->posc_j(P.Y[n]);
        k=p->posc_k(P.Z[n]);
        }
        
        cellSum(i,j,k) += P.ParcelFactor;
        bedch(i,j) += P.ParcelFactor;
    }
    
    pgc->gcsl_start4(p,bedch,1);
}

void partres::cellSum_full_update(lexer *p, ghostcell *pgc, sediment_fdm *s, int mode)
{
    ALOOP
    cellSum(i,j,k) = 0.0;
    pgc->start4a(p,cellSum,1);
    
    for(size_t n=0;n<P.index;n++)
    if(P.Flag[n]==ACTIVE)
    {
        if(mode==1)
        {
        i=p->posc_i(P.XRK1[n]);
        j=p->posc_j(P.YRK1[n]);
        k=p->posc_k(P.ZRK1[n]);
        
        Sx = (p->XN[IP1] - P.XRK1[n])/(p->XN[IP1] - p->XN[IP]);
        Sy = (p->YN[JP1] - P.YRK1[n])/(p->YN[JP1] - p->YN[JP]);
        Sz = (p->ZN[KP1] - P.ZRK1[n])/(p->ZN[KP1] - p->ZN[KP]);
        }
        
        if(mode==2)
        {
        i=p->posc_i(P.X[n]);
        j=p->posc_j(P.Y[n]);
        k=p->posc_k(P.Z[n]);
        
        Sx = (p->XN[IP1] - P.X[n])/(p->XN[IP1] - p->XN[IP]);
        Sy = (p->YN[JP1] - P.Y[n])/(p->YN[JP1] - p->YN[JP]);
        Sz = (p->ZN[KP1] - P.Z[n])/(p->ZN[KP1] - p->ZN[KP]);
        }
        
        cellSum(i,j,k) += P.ParcelFactor * Sx*Sy*Sz;
        cellSum(i+1,j,k) += P.ParcelFactor * (1.0-Sx)*Sy*Sz;
        cellSum(i+1,j+1,k) += P.ParcelFactor * (1.0-Sx)*(1.0-Sy)*Sz;
        cellSum(i,j+1,k) += P.ParcelFactor * Sx*(1.0-Sy)*Sz;
        cellSum(i,j,k+1) += P.ParcelFactor * Sx*Sy*(1.0-Sz);
        cellSum(i+1,j,k+1) += P.ParcelFactor * (1.0-Sx)*Sy*(1.0-Sz);
        cellSum(i+1,j+1,k+1) += P.ParcelFactor * (1.0-Sx)*(1.0-Sy)*(1.0-Sz);
        cellSum(i,j+1,k+1) += P.ParcelFactor * Sx*(1.0-Sy)*(1.0-Sz);
    }
    
    pgc->start4a_sum(p,cellSum,1);
}




