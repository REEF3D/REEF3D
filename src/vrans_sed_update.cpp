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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"vrans_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void vrans_f::sed_update(lexer *p, fdm *a, ghostcell *pgc)
{
    ALOOP
		if(a->topo(i,j,k)<0.0) // Topo
		{
			a->porosity(i,j,k) = p->S24; //porosity
			porpart(i,j,k) = p->S20;  //d50
			alpha(i,j,k) = p->S26_a;  //alpha
			beta(i,j,k) = p->S26_b;    //beta
		}
		else // Not topo
		{
			a->porosity(i,j,k)=1.0;
			porpart(i,j,k)=0.01;
			alpha(i,j,k)=0.0;
			beta(i,j,k)=0.0;
		}
    
    
    pgc->start4a(p,a->porosity,1);
	pgc->start4a(p,porpart,1);
	pgc->start4a(p,alpha,1);
	pgc->start4a(p,beta,1);
}

