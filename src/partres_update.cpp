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
#include"fdm.h"
#include"ghostcell.h"
#include"vrans.h"
#include"sediment_fdm.h"

void partres::update(lexer *p, fdm *a, ghostcell *pgc, sediment_fdm *s, field &por, field &d50)
{
    
    ALOOP
	{
	por(i,j,k)=p->S24;
	d50(i,j,k)=p->S20;
    
    //a->visc(i,j,k) = p->W2*pow((1.0 + 0.5*p->W2*Ts(i,j,k))/(1.0 - Ts(i,j,k)/0.6) ,2.0);
	}
    
    
    pgc->start4a(p,por,1);
    pgc->start4a(p,d50,1);
    
}