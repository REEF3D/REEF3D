/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"particle_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<math.h>

void particle_f::remove(lexer* p, fdm* a, ghostcell* pgc)
{
    bool inBounds=false;
    removed=0;

    PARTLOOP
        if(PP.Flag[n]>0)
        {
            i = p->posc_i(PP.X[n]);
            j = p->posc_j(PP.Y[n]);
            k = p->posc_k(PP.Z[n]);

            inBounds=minboundcheck(p,i,j,k,1);
            if (inBounds)
                inBounds=maxboundcheck(p,i,j,k,1);

			// remove out of bounds particles
            if(!inBounds)
            {
                PP.erase(n);
                removed++;
            }
        }
}
