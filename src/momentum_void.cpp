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

#include"momentum_void.h"

momentum_void::momentum_void()
{
}

momentum_void::~momentum_void()
{
}

void momentum_void::start(lexer *p, fdm* a, ghostcell* pgc, vrans* pvrans, sixdof *p6dof, vector<net*>& pnet)
{
}

void momentum_void::utimesave(lexer *p, fdm *a, ghostcell *pgc)
{
}

void momentum_void::vtimesave(lexer *p, fdm *a, ghostcell *pgc)
{
}

void momentum_void::wtimesave(lexer *p, fdm *a, ghostcell *pgc)
{
}

