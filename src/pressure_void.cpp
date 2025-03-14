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

#include"pressure_void.h"

void pressure_void::start(lexer*, fdm*, ghostcell*, ioflow*, solver*, field&, field&, field&, double)
{
}

void pressure_void::ucorr(lexer*, fdm*, field&, double)
{    
}

void pressure_void::vcorr(lexer*, fdm*, field&, double)
{     
}

void pressure_void::wcorr(lexer*, fdm*, field&, double)
{    
}

void pressure_void::upgrad(lexer*, fdm*, slice&, slice&)
{
}

void pressure_void::vpgrad(lexer*, fdm*, slice&, slice&)
{
}

void pressure_void::wpgrad(lexer*, fdm*, slice&, slice&)
{
}

void pressure_void::ini(lexer*, fdm*, ghostcell*)
{
}
