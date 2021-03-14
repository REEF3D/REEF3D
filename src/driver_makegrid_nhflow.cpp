/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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
--------------------------------------------------------------------*/

#include"driver.h"
#include"lexer.h"
#include"ghostcell.h"

void driver::makegrid_nhflow(lexer *p, ghostcell *pgc)
{	
    int q;
    
// flag7
    p->Iarray(p->flag7,p->imax*p->jmax*(p->kmax+2));
    
    for(i=0;i<p->imax*p->jmax*(p->kmax+2);++i)
    p->flag7[i]=-10;
    
    BASELOOP
    {
        p->flag7[FIJK]=p->flag4[IJK];
    }
    
    k=p->knoz;
    SLICEBASELOOP
    p->flag7[FIJK] = p->flag7[FIJKm1];
    
    
}
