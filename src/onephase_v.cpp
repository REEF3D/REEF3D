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
along with this program; if not, see <http://www.gnu.org/liceonephases/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"onephase_v.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ioflow.h"


onephase_v::onephase_v(lexer *p, fdm *a, ghostcell *pgc) 
{
}

onephase_v::~onephase_v()
{
}

void onephase_v::update(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow)
{    
}

void onephase_v::ini(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow)
{
}

void onephase_v::uvel(lexer*, fdm*, ghostcell*, field&)
{
}

void onephase_v::vvel(lexer*, fdm*, ghostcell*, field&)
{
}

void onephase_v::wvel(lexer*, fdm*, ghostcell*, field&)
{
}
