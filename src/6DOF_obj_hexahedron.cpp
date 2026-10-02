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

void sixdof_obj::hexahedron(lexer *p, ghostcell *pgc, int id)
{
    const double v[24] = {p->X164_x1[id], p->X164_y1[id], p->X164_z1[id], p->X164_x2[id], p->X164_y2[id], p->X164_z2[id], p->X164_x3[id], p->X164_y3[id], p->X164_z3[id], p->X164_x4[id], p->X164_y4[id], p->X164_z4[id], p->X164_x5[id], p->X164_y5[id], p->X164_z5[id], p->X164_x6[id], p->X164_y6[id], p->X164_z6[id], p->X164_x7[id], p->X164_y7[id], p->X164_z7[id], p->X164_x8[id], p->X164_y8[id], p->X164_z8[id]};
    
	tstart[entity_count]=tricount;
    
    geo_primitive::hexahedron(tri_x,tri_y,tri_z,tricount,v);
    
    tend[entity_count]=tricount;
}
