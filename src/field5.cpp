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

#include"field5.h"
#include"lexer.h"

field5::field5(lexer *p)
{
    imin=p->imin;
    imax=p->imax;
    jmin=p->jmin;
    jmax=p->jmax;
    kmin=p->kmin;
    kmax=p->kmax;

	fieldalloc(p);
}

field5::~field5()
{
    delete [] V;
}

void field5::fieldalloc(lexer* p)
{
	int gridsize = imax*jmax*kmax;
	p->Darray(V,gridsize);
}

void field5::resize(lexer* p)
{
}

void field5::dealloc(lexer* p)
{
	delete [] V;
}

double & field5::operator[](int n)
{
	return V[n];
}

double & field5::operator()(int ii, int jj, int kk)
{
    iter=(ii-imin)*jmax*kmax + (jj-jmin)*kmax + kk-kmin;
	return V[iter];

}


