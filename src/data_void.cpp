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

#include"data_void.h"

data_void::data_void()
{
}

void data_void::start(lexer*, fdm*, ghostcell*)
{
}

void data_void::print_3D(lexer*, fdm*, ghostcell*, ofstream&)
{
}

void data_void::name_pvtu(lexer*, fdm*, ghostcell*, ofstream&)
{
}

void data_void::name_vtu(lexer*, fdm*, ghostcell*, ofstream&, int*, int&)
{
}

void data_void::offset_vtu(lexer*, fdm*, ghostcell*, ofstream&, int*, int&)
{
}
