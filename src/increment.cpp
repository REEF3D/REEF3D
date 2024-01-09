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

#include"increment.h"
#include"fdm.h"
#include"fdm2D.h"

increment::increment()
{
	pip=0;
    marge=5;
}

increment::~increment()
{
}

int increment::i,increment::j,increment::k,increment::n,increment::h,increment::innercounter,increment::pip;
int increment::marge;
fdm* increment::aa;
fdm2D* increment::bb;

//,increment::l
