/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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
--------------------------------------------------------------------*/

#include"patch_obj.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

patch_obj::patch_obj(lexer *p, int ID_ini) 
{
    ID = ID_ini;
    
    gcb_count=0;
    
}

patch_obj::~patch_obj()
{
}

void patch_obj::patch_obj_ini(lexer *p, ghostcell *pgc)
{
}

void patch_obj::patch_obj_gcb_generate(lexer *p, ghostcell *pgc)
{
    p->Iarray(gcb,gcb_count,3);
}