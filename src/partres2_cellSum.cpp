/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"partres2.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment_fdm.h"

void partres2::cellSum_update(lexer *p, ghostcell *pgc, sediment_fdm *s, int mode)
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

void partres2::cellSum_full_update(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
    BLOOP
    cellSum(i,j,k) = 0;
    
    for(size_t n=0;n<P.index;n++)
    if(P.Flag[n]==ACTIVE)
    {
        i=p->posc_i(P.X[n]);
        j=p->posc_j(P.Y[n]);
        k=p->posc_k(P.Z[n]);
        
        cellSum(i,j,k) += P.ParcelFactor;
    }
    
}