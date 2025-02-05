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

#include"part.h"
#include"lexer.h"
#include"ghostcell.h"
#include"slice.h"

void part::xchange_fillback(lexer *p, ghostcell *pgc, double *F)
{
    index_empty = index_empty0;
    
    // fill recv into F
    for(n=0;n<6;++n)
    for(q=0;q<recvnum[n];++q)
    {
    F[Empty[index_empty]] = recv[n][q]; 
    --index_empty;
    }
}

void part::xchange_fillback_flag(lexer *p, ghostcell *pgc, slice &bedch, int mode)
{
    index_empty = index_empty0;
    
    // fill recv into F
    n=0;
    for(int qn=0;qn<6;++qn)
    for(q=0;q<recvnum[qn];++q)
    {
    n=Empty[index_empty];
    
    // flag
    Flag[n] = ACTIVE; 
    
    // bedch
    if(mode==1)
    {
    i=p->posc_i(XRK1[n]);
    j=p->posc_j(YRK1[n]);
    }
            
    if(mode==2)
    {
    i=p->posc_i(X[n]);
    j=p->posc_j(Y[n]);
    }
            
    bedch(i,j) += ParcelFactor;
    
    --index_empty;
    ++n;
    }
}



