/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"cpt.h"
#include"lexer.h"

cpt::cpt()
{
}

cpt::~cpt()
{
    delete [] n;
	delete [] s;
	delete [] w;
	delete [] e;
	delete [] t;
	delete [] b;
	delete [] p;
}

void cpt::allocate(lexer *pp)
{
    pp->Iarray(n,pp->veclength);
	pp->Iarray(s,pp->veclength);
	pp->Iarray(e,pp->veclength);
	pp->Iarray(w,pp->veclength);
	pp->Iarray(t,pp->veclength);
	pp->Iarray(b,pp->veclength);
	pp->Iarray(p,pp->veclength);
}

void cpt::resize(lexer *pp, int size_old, int size_new)
{
    pp->Iresize(n,size_old,size_new);
    pp->Iresize(s,size_old,size_new);
    pp->Iresize(e,size_old,size_new);
    pp->Iresize(w,size_old,size_new);
    pp->Iresize(t,size_old,size_new);
    pp->Iresize(b,size_old,size_new);
    pp->Iresize(p,size_old,size_new);
}

