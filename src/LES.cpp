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

#include"LES.h"
#include"lexer.h"
#include"fdm.h"

LES::LES(lexer* p, fdm* a) : les_io(p,a)
{
}

LES::~LES()
{
}

void LES::isource(lexer* p, fdm* a)
{
	ULOOP
	a->F(i,j,k)=0.0;
}

void LES::jsource(lexer* p, fdm* a)
{
	VLOOP
	a->G(i,j,k)=0.0;
}

void LES::ksource(lexer* p, fdm* a)
{
	WLOOP
	a->H(i,j,k)=0.0;
}







