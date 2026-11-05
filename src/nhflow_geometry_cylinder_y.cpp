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

#include"nhflow_geometry.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"geo_primitive.h"

void nhflow_geometry::cylinder_y(lexer *p, ghostcell *pgc, int id)
{
    const int snum = geo_primitive::segments(cyly_r[id],0.5*(DSM),40);
    
	tstart[entity_count]=tricount;
    
    geo_primitive::cylinder_y(tri_x,tri_y,tri_z,tricount,cyly_xc[id],cyly_zc[id],cyly_ys[id],cyly_ye[id],cyly_r[id],snum);
	
	tend[entity_count]=tricount;
}
