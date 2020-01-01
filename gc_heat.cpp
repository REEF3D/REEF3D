/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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

#include"lexer.h"
#include"ghostcell.h"
#include"field.h"
#include"vec.h"
#include"fdm.h"

void ghostcell::heatbc(lexer *p, field& f, int gcv, int bc, int cs)
{
	if(cs==1)
	for(q=0;q<=margin;++q)
	f(i-q,j,k)=p->H61_T;

	if(cs==2)
	for(q=0;q<=margin;++q)
	f(i,j+q,k)=p->H62_T;

	if(cs==3)
	for(q=0;q<=margin;++q)
	f(i,j-q,k)=p->H63_T;

	if(cs==4)
	for(q=0;q<=margin;++q)
	f(i+q,j,k)=p->H64_T;

	if(cs==5)
	for(q=0;q<=margin;++q)
	f(i,j,k-q)=p->H65_T;

	if(cs==6)
	for(q=0;q<=margin;++q)
	f(i,j,k+q)=p->H66_T;
}
