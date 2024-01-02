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

#include"vrans_veg.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

vrans_veg::vrans_veg(lexer *p, ghostcell *pgc) : Cval(p->B264), N(p), D(p), Cd(p), un(p), vn(p), wn(p)
{
	//initialize(p,a,pgc);
}

vrans_veg::~vrans_veg()
{
}

void vrans_veg::veltimesave(lexer *p, fdm *a, ghostcell *pgc)
{
    ULOOP
    un(i,j,k) = a->u(i,j,k);
    
    VLOOP
    vn(i,j,k) = a->v(i,j,k);
    
    WLOOP
    wn(i,j,k) = a->w(i,j,k);
}

void vrans_veg::sed_update(lexer *p, fdm *a, ghostcell *pgc)
{
}
