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

#include"ghostcell.h"
#include"lexer.h"
#include"field.h"

// diagonal ghost cells (Y 40): dgc[n] = {i,j,k, di,dj,dk, bc}, bc 1 copy of (i,j,k), bc 2 zero
void ghostcell::dgcpol(lexer* p, field& f, int **dgc, int dgc_count)
{
    int di,dj,dk,bc;
    
    for(n=0;n<dgc_count;++n)
    {
        i=dgc[n][0];
        j=dgc[n][1];
        k=dgc[n][2];
        
        di=dgc[n][3];
        dj=dgc[n][4];
        dk=dgc[n][5];
        
        bc=dgc[n][6];
        
        if(bc==1)
        f(i+di,j+dj,k+dk) = f(i,j,k);  

        if(bc==2)
        f(i+di,j+dj,k+dk) = 0.0;        
    }
}
