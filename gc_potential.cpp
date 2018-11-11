/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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

#include"ghostcell.h"
#include"field.h"
#include"lexer.h"
#include"fdm.h"

void ghostcell::potentialbc(lexer *p, field& f, int bc, int cs)
{
    if(cs==1)
	f(i-1,j,k) =  0.0;//-a->u(i-1,j,k)*deltax + f(i,j,k);

	if(cs==4)
	f(i+1,j,k) =  a->u(i,j,k)*deltax + f(i,j,k);
}

