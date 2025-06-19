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

#include"sliceint4.h"
#include"lexer.h"
#include"fdm.h"

sliceint4::sliceint4(lexer *p)
{
    imin=p->imin;
    imax=p->imax;
    jmin=p->jmin;
    jmax=p->jmax;
    
	fieldalloc(p);
    
	pp=p;
}

sliceint4::~sliceint4()
{
	delete [ ] V;
}

void sliceint4::fieldalloc(lexer* p)
{
	int gridsize = imax*jmax;
	p->Iarray(V,gridsize);
}

void sliceint4::dealloc(lexer* p)
{
	delete [ ] V;
}

void sliceint4::resize(lexer* p)
{
}

int & sliceint4::operator()(int ii, int jj)
{		
	return V[(ii-imin)*jmax + (jj-jmin)];
}

