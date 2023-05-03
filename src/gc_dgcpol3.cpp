/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

void ghostcell::dgcpol3(lexer* p,field& f,int gcv)
{
    int di,dj,dk,bc;
    
    for(n=0;n<p->dgc3_count;++n)
    {
        i=p->dgc3[n][0];
        j=p->dgc3[n][1];
        k=p->dgc3[n][2];
        
        
        di=p->dgc3[n][3];
        dj=p->dgc3[n][4];
        dk=p->dgc3[n][5];
        
        bc=p->dgc3[n][6];
        
        if(bc==1)
        f(i+di,j+dj,k+dk) = f(i,j,k);  

        if(bc==2)
        f(i+di,j+dj,k+dk) = 0.0;        
        
    }
}
