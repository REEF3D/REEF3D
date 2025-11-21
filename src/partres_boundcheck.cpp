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

void partres::boundcheck(lexer *p, int mode)
{
    bool inBounds;

    for(n=0;n<P.index;++n)
    if(P.Flag[n]==ACTIVE)
    {
        if(mode==1)
        {
            i=p->posc_i(P.XRK1[n]);
            j=p->posc_j(P.YRK1[n]);
            k=p->posc_k(P.ZRK1[n]);
        }
        else if(mode==2)
        {
            i=p->posc_i(P.X[n]);
            j=p->posc_j(P.Y[n]);
            k=p->posc_k(P.Z[n]);
        }

        inBounds=boundaries.globalminboundcheck(p,i,j,k);
        if(inBounds)
            inBounds=boundaries.globalmaxboundcheck(p,i,j,k);

        // remove out of bounds particles
        if(inBounds==0)
            P.remove(n);
    }
}
