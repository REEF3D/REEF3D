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

void sixdof_obj::cylinder_x(lexer *p, ghostcell *pgc, int id)
{
    // segment length 0.75 dx (was 0.75*U*dx, i.e. 1/(0.75 dx) segments for any radius)
	const int snum = geo_primitive::segments(p->X131_rad,0.75*p->dx,8);
	
	tstart[entity_count]=tricount;
	
	geo_primitive::cylinder_x(tri_x,tri_y,tri_z,tricount,p->X131_yc,p->X131_zc,p->X131_xc-0.5*p->X131_h,p->X131_xc+0.5*p->X131_h,p->X131_rad,snum);
	
	tend[entity_count]=tricount;
}

void sixdof_obj::cylinder_y(lexer *p, ghostcell *pgc, int id)
{
    // segment length 0.75 dx (was 0.75*U*dx, i.e. 1/(0.75 dx) segments for any radius)
	const int snum = geo_primitive::segments(p->X132_rad,0.75*p->dx,8);
	
	tstart[entity_count]=tricount;
	
	geo_primitive::cylinder_y(tri_x,tri_y,tri_z,tricount,p->X132_xc,p->X132_zc,p->X132_yc-0.5*p->X132_h,p->X132_yc+0.5*p->X132_h,p->X132_rad,snum);
	
	tend[entity_count]=tricount;
}

void sixdof_obj::cylinder_z(lexer *p, ghostcell *pgc, int id)
{
    // segment length 0.75 dx (was 0.75*U*dx, i.e. 1/(0.75 dx) segments for any radius)
	const int snum = geo_primitive::segments(p->X133_rad,0.75*p->dx,8);
	
	tstart[entity_count]=tricount;
	
	geo_primitive::cylinder_z(tri_x,tri_y,tri_z,tricount,p->X133_xc,p->X133_yc,p->X133_zc-0.5*p->X133_h,p->X133_zc+0.5*p->X133_h,p->X133_rad,snum);
	
	tend[entity_count]=tricount;
}
