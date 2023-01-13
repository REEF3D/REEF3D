/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"patchBC_void.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"patch_obj.h"

patchBC_void::patchBC_void(lexer *p) 
{
    obj_count=0;
}

patchBC_void::~patchBC_void()
{
}

void patchBC_void::patchBC_ini(lexer *p, ghostcell *pgc)
{
} 

void patchBC_void::patchBC_ioflow(lexer *p, fdm *a, ghostcell *pgc, field &u, field &v, field &w)
{
} 

void patchBC_void::patchBC_rkioflow(lexer *p, fdm *a, ghostcell *pgc, field &u, field &v, field &w)
{
} 

void patchBC_void::patchBC_discharge(lexer *p, fdm *a, ghostcell *pgc)
{
} 

void patchBC_void::patchBC_pressure(lexer *p, fdm *a, ghostcell *pgc, field &press)
{
} 

void patchBC_void::patchBC_waterlevel(lexer *p, fdm *a, ghostcell *pgc, field &phi)
{
} 

void patchBC_void::patchBC_ioflow2D(lexer *p, ghostcell*, slice&, slice&, slice&, slice&)
{
}
void patchBC_void::patchBC_rkioflow2D(lexer *p, ghostcell*, slice&, slice&, slice&, slice&)
{
}

void patchBC_void::patchBC_discharge2D(lexer *p, fdm2D*, ghostcell*, slice&, slice&, slice&, slice&)
{
}

void patchBC_void::patchBC_pressure2D(lexer*, ghostcell*, slice&)
{
}

void patchBC_void::patchBC_pressure2D_ugrad(lexer*, fdm2D*, slice&,slice&)
{
}

void patchBC_void::patchBC_pressure2D_vgrad(lexer*, fdm2D*, slice&, slice&)
{
}

void patchBC_void::patchBC_waterlevel2D(lexer*, fdm2D*, ghostcell*, slice&)
{
}

void patchBC_void::patchBC_loop2D(lexer*, fdm2D*, int&, int&, int&, int&)
{
}
