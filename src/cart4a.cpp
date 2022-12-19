/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2022 Hans Bihs

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

#include"cart4a.h"

using namespace std;

cart4a::cart4a(lexer *p)
{
}

cart4a::~cart4a()
{
	delete [] ggc;
	delete [] ggcmem;
}

//-----------------------------------------------------------------------------

int** cart4a::ggc;
int* cart4a::ggcmem;
int cart4a::ggccount;
int cart4a::ggcsize;

