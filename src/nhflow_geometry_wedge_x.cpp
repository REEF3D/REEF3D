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

void nhflow_geometry::wedge_x(lexer *p, ghostcell *pgc, int id)
{
    xs = wedgex_xs[id];
    xe = wedgex_xe[id];
	
    ys = wedgex_ys[id];
    ye = wedgex_ye[id];

    zs = wedgex_zs[id];
    ze = wedgex_ze[id];
	
	tstart[entity_count]=tricount;
    
    geo_primitive::wedge_x(tri_x,tri_y,tri_z,tricount,xs,xe,ys,ye,zs,ze);
	
	tend[entity_count]=tricount;
}
