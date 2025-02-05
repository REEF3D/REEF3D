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

void ghostcell::solid_forcing_lsm(lexer *p, fdm *a, field &f)
{

    if(p->X48==1)
    GCDF4LOOP
    {
    i=p->gcdf4[n][0];
    j=p->gcdf4[n][1];
    k=p->gcdf4[n][2];

        if(p->gcdf4[n][3]==1)
        {
        f(i-1,j,k)=f(i,j,k);
        f(i-2,j,k)=f(i,j,k);
        }

        if(p->gcdf4[n][3]==4)
        {
        f(i+1,j,k)=f(i,j,k);
        f(i+2,j,k)=f(i,j,k);
        }

        if(p->gcdf4[n][3]==3)
        {
        f(i,j-1,k)=f(i,j,k);
        f(i,j-2,k)=f(i,j,k);
        }

        if(p->gcdf4[n][3]==2)
        {
        f(i,j+1,k)=f(i,j,k);
        f(i,j+2,k)=f(i,j,k);
        }

        if(p->gcdf4[n][3]==5)
        {
        f(i,j,k-1)=f(i,j,k);
        f(i,j,k-2)=f(i,j,k);
        }

        if(p->gcdf4[n][3]==6)
        {
        f(i,j,k+1)=f(i,j,k);
        f(i,j,k+2)=f(i,j,k);
        }
	}
}
