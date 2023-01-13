/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"convection_void.h"
#include"lexer.h"
#include"fdm.h"

convection_void::convection_void (lexer *p)
{
}

convection_void::~convection_void()
{
}

void convection_void::start(lexer* p, fdm* a, field& b, int ipol, field& uvel, field& vvel, field& wvel)
{
	int count=0;
	LOOP
	{
	 a->M.p[count] = 0.0;
	 
	 a->M.s[count] = 0.0;
	 a->M.n[count] = 0.0;
	 
	 a->M.e[count] = 0.0;
	 a->M.w[count] = 0.0;
	 
	 a->M.b[count] = 0.0;
	 a->M.t[count] = 0.0;
	 
	 ++count;
	}
}
