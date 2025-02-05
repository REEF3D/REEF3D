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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"field2.h"
#include"lexer.h"
#include"fdm.h"

field2::field2(lexer *p)
{
    imin=p->imin;
    imax=p->imax;
    jmin=p->jmin;
    jmax=p->jmax;
    kmin=p->kmin;
    kmax=p->kmax;

	fieldalloc(p);

	pp=p;
}

field2::~field2()
{
	delete [ ] V;
}

void field2::fieldalloc(lexer* p)
{
	int gridsize = p->imax*p->jmax*p->kmax;
	p->Darray(V,gridsize);
}

void field2::dealloc(lexer* p)
{
	delete [ ] V;
}

void field2::resize(lexer* p)
{
}

double & field2::operator[](int n)
{
	return V[n];
}

double & field2::operator()(int ii, int jj, int kk)
{	
	return V[(ii-imin)*jmax*kmax + (jj-jmin)*kmax + kk-kmin];
}
