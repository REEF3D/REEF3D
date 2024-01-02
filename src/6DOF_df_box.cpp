/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"6DOF_df_object.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void sixdof_df_object::box(lexer *p, fdm *a, ghostcell *pgc, int id)
{
	xs = p->X110_xs[id];
    xe = p->X110_xe[id];
	
    ys = p->X110_ys[id];
    ye = p->X110_ye[id];

    zs = p->X110_zs[id];
    ze = p->X110_ze[id];    
	
	// Face 3
	// Tri 1
	tstart[entity_count]=tricount;
	
	tri_x[tricount][0] = xs;
	tri_x[tricount][1] = xe;
	tri_x[tricount][2] = xe;

	tri_y[tricount][0] = ys;
	tri_y[tricount][1] = ys;
	tri_y[tricount][2] = ys;

	tri_z[tricount][0] = zs;
	tri_z[tricount][1] = zs;
	tri_z[tricount][2] = ze;
	++tricount;

	// Tri 2
	tri_x[tricount][0] = xe;
	tri_x[tricount][1] = xs;
	tri_x[tricount][2] = xs;

	tri_y[tricount][0] = ys;
	tri_y[tricount][1] = ys;
	tri_y[tricount][2] = ys;

	tri_z[tricount][0] = ze;
	tri_z[tricount][1] = ze;
	tri_z[tricount][2] = zs;
	++tricount;

	// Face 4
	// Tri 3
	tri_x[tricount][0] = xe;
	tri_x[tricount][1] = xe;
	tri_x[tricount][2] = xe;

	tri_y[tricount][0] = ye;
	tri_y[tricount][1] = ys;
	tri_y[tricount][2] = ys;

	tri_z[tricount][0] = ze;
	tri_z[tricount][1] = ze;
	tri_z[tricount][2] = zs;
	++tricount;

	// Tri 4
	tri_x[tricount][0] = xe;
	tri_x[tricount][1] = xe;
	tri_x[tricount][2] = xe;

	tri_y[tricount][0] = ye;
	tri_y[tricount][1] = ye;
	tri_y[tricount][2] = ys;

	tri_z[tricount][0] = zs;
	tri_z[tricount][1] = ze;
	tri_z[tricount][2] = zs;
	++tricount;

	// Face 1
	// Tri 5
	tri_x[tricount][0] = xs;
	tri_x[tricount][1] = xs;
	tri_x[tricount][2] = xs;

	tri_y[tricount][0] = ye;
	tri_y[tricount][1] = ye;
	tri_y[tricount][2] = ys;

	tri_z[tricount][0] = ze;
	tri_z[tricount][1] = zs;
	tri_z[tricount][2] = zs;
	++tricount;

	// Tri 6
	tri_x[tricount][0] = xs;
	tri_x[tricount][1] = xs;
	tri_x[tricount][2] = xs;

	tri_y[tricount][0] = ys;
	tri_y[tricount][1] = ye;
	tri_y[tricount][2] = ys;

	tri_z[tricount][0] = ze;
	tri_z[tricount][1] = ze;
	tri_z[tricount][2] = zs;
	++tricount;
	
	// Face 2
	// Tri 7
	tri_x[tricount][0] = xe;
	tri_x[tricount][1] = xe;
	tri_x[tricount][2] = xs;

	tri_y[tricount][0] = ye;
	tri_y[tricount][1] = ye;
	tri_y[tricount][2] = ye;

	tri_z[tricount][0] = ze;
	tri_z[tricount][1] = zs;
	tri_z[tricount][2] = zs;
	++tricount;

	// Tri 8
	tri_x[tricount][0] = xs;
	tri_x[tricount][1] = xe;
	tri_x[tricount][2] = xs;

	tri_y[tricount][0] = ye;
	tri_y[tricount][1] = ye;
	tri_y[tricount][2] = ye;

	tri_z[tricount][0] = ze;
	tri_z[tricount][1] = ze;
	tri_z[tricount][2] = zs;
	++tricount;

	// Face 5
	// Tri 9
	tri_x[tricount][0] = xe;
	tri_x[tricount][1] = xe;
	tri_x[tricount][2] = xs;

	tri_y[tricount][0] = ye;
	tri_y[tricount][1] = ys;
	tri_y[tricount][2] = ys;

	tri_z[tricount][0] = zs;
	tri_z[tricount][1] = zs;
	tri_z[tricount][2] = zs;
	++tricount;

	// Tri 10
	tri_x[tricount][0] = xs;
	tri_x[tricount][1] = xe;
	tri_x[tricount][2] = xs;

	tri_y[tricount][0] = ye;
	tri_y[tricount][1] = ye;
	tri_y[tricount][2] = ys;

	tri_z[tricount][0] = zs;
	tri_z[tricount][1] = zs;
	tri_z[tricount][2] = zs;
	++tricount;

	// Face 6
	// Tri 11
	tri_x[tricount][0] = xs;
	tri_x[tricount][1] = xe;
	tri_x[tricount][2] = xe;

	tri_y[tricount][0] = ys;
	tri_y[tricount][1] = ys;
	tri_y[tricount][2] = ye;

	tri_z[tricount][0] = ze;
	tri_z[tricount][1] = ze;
	tri_z[tricount][2] = ze;
	++tricount;

	// Tri 12
	tri_x[tricount][0] = xs;
	tri_x[tricount][1] = xe;
	tri_x[tricount][2] = xs;

	tri_y[tricount][0] = ys;
	tri_y[tricount][1] = ye;
	tri_y[tricount][2] = ye;

	tri_z[tricount][0] = ze;
	tri_z[tricount][1] = ze;
	tri_z[tricount][2] = ze;
	++tricount;
	 
	tend[entity_count]=tricount;
}

