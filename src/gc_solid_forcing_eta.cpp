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
Authors: Hans Bihs
--------------------------------------------------------------------*/

#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"

void ghostcell::solid_forcing_eta(lexer *p, slice &f)
{
    
    //if(p->X48==1)
    GCSLDF4LOOP
    {
    i=p->gcsldf4[n][0];
    j=p->gcsldf4[n][1];

        if(p->gcsldf4[n][3]==1)
        {
        f(i-1,j)=f(i,j);
        f(i-2,j)=f(i,j);
        }

        if(p->gcsldf4[n][3]==4)
        {
        f(i+1,j)=f(i,j);
        f(i+2,j)=f(i,j);
        }

        if(p->gcsldf4[n][3]==3)
        {
        f(i,j-1)=f(i,j);
        f(i,j-2)=f(i,j);
        }

        if(p->gcsldf4[n][3]==2)
        {
        f(i,j+1)=f(i,j);
        f(i,j+2)=f(i,j);
        }

	}
}
