/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#include"6DOF_obj.h"
#include"lexer.h"
#include"ghostcell.h"
#include"geo_primitive.h"

void sixdof_obj::wedge_sym(lexer *p, ghostcell *pgc, int id)
{
	xs = p->X153_xs;
    xe = p->X153_xe;
	
    ys = p->X153_ys;
    ye = p->X153_ye;

    zs = p->X153_zs;
    ze = p->X153_ze;  
    
	tstart[entity_count]=tricount;
    
    geo_primitive::wedge_sym(tri_x,tri_y,tri_z,tricount,xs,xe,ys,ye,zs,ze);
    
	tend[entity_count]=tricount;
}
