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

void particle_f::parcount(lexer* p, fdm* a, ghostcell* pgc)
{		
	LOOP
	posnum(i,j,k)=0.0;

	
	pgc->start4(p,posnum,1);
		
	// POS
    for(n=0;n<posactive;++n)
    if(posflag[n]>0)
    {
        i = p->posc_i(pos[n][0]);
        j = p->posc_j(pos[n][1]);
        k = p->posc_k(pos[n][2]);

    posnum(i,j,k)+=1.0;
    }


}
