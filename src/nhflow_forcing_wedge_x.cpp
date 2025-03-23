/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"nhflow_forcing.h"
#include"lexer.h"
#include"ghostcell.h"

void nhflow_forcing::wedge_x(lexer *p, ghostcell *pgc, int id)
{
    xs = p->A587_xs[id];
    xe = p->A587_xe[id];
	
    ys = p->A587_ys[id];
    ye = p->A587_ye[id];

    zs = p->A587_zs[id];
    ze = p->A587_ze[id];
    
	// Face 3
	// Tri 1
	tstart[entity_count]=tricount;
	
	if(zs<ze)
	{
	// sides
	tri_x[tricount][0] = xs;
	tri_y[tricount][0] = ys;
	tri_z[tricount][0] = zs;

	tri_x[tricount][1] = xe;
	tri_y[tricount][1] = ys;
	tri_z[tricount][1] = zs;

	tri_x[tricount][2] = xe;
	tri_y[tricount][2] = ys;
	tri_z[tricount][2] = ze;
	++tricount;

	tri_x[tricount][0] = xs;
	tri_y[tricount][0] = ye;
	tri_z[tricount][0] = zs;

	tri_x[tricount][1] = xe;
	tri_y[tricount][1] = ye;
	tri_z[tricount][1] = zs;

	tri_x[tricount][2] = xe;
	tri_y[tricount][2] = ye;
	tri_z[tricount][2] = ze;
	++tricount;

	// bottom
	tri_x[tricount][0] = xs;
	tri_y[tricount][0] = ys;
	tri_z[tricount][0] = zs;

	tri_x[tricount][1] = xe;
	tri_y[tricount][1] = ys;
	tri_z[tricount][1] = zs;

	tri_x[tricount][2] = xe;
	tri_y[tricount][2] = ye;
	tri_z[tricount][2] = zs;
	++tricount;

	tri_x[tricount][0] = xs;
	tri_y[tricount][0] = ys;
	tri_z[tricount][0] = zs;

	tri_x[tricount][1] = xe;
	tri_y[tricount][1] = ye;
	tri_z[tricount][1] = zs;

	tri_x[tricount][2] = xs;
	tri_y[tricount][2] = ye;
	tri_z[tricount][2] = zs;
	++tricount;

	// front
	tri_x[tricount][0] = xe;
	tri_y[tricount][0] = ys;
	tri_z[tricount][0] = zs;

	tri_x[tricount][1] = xe;
	tri_y[tricount][1] = ye;
	tri_z[tricount][1] = zs;

	tri_x[tricount][2] = xe;
	tri_y[tricount][2] = ys;
	tri_z[tricount][2] = ze;
	++tricount;

	tri_x[tricount][0] = xe;
	tri_y[tricount][0] = ye;
	tri_z[tricount][0] = zs;

	tri_x[tricount][1] = xe;
	tri_y[tricount][1] = ye;
	tri_z[tricount][1] = ze;

	tri_x[tricount][2] = xe;
	tri_y[tricount][2] = ys;
	tri_z[tricount][2] = ze;
	++tricount;

	// top
	tri_x[tricount][0] = xs;
	tri_y[tricount][0] = ys;
	tri_z[tricount][0] = zs;

	tri_x[tricount][1] = xe;
	tri_y[tricount][1] = ys;
	tri_z[tricount][1] = ze;

	tri_x[tricount][2] = xe;
	tri_y[tricount][2] = ye;
	tri_z[tricount][2] = ze;
	++tricount;

	tri_x[tricount][0] = xs;
	tri_y[tricount][0] = ys;
	tri_z[tricount][0] = zs;

	tri_x[tricount][1] = xe;
	tri_y[tricount][1] = ye;
	tri_z[tricount][1] = ze;

	tri_x[tricount][2] = xs;
	tri_y[tricount][2] = ye;
	tri_z[tricount][2] = zs;
	++tricount;


	}

	if(zs>=ze)
	{
	// sides
	tri_x[tricount][0] = xs;
	tri_y[tricount][0] = ys;
	tri_z[tricount][0] = ze;

	tri_x[tricount][1] = xe;
	tri_y[tricount][1] = ys;
	tri_z[tricount][1] = ze;

	tri_x[tricount][2] = xs;
	tri_y[tricount][2] = ys;
	tri_z[tricount][2] = zs;
	++tricount;

	tri_x[tricount][0] = xs;
	tri_y[tricount][0] = ye;
	tri_z[tricount][0] = ze;

	tri_x[tricount][1] = xe;
	tri_y[tricount][1] = ye;
	tri_z[tricount][1] = ze;

	tri_x[tricount][2] = xs;
	tri_y[tricount][2] = ye;
	tri_z[tricount][2] = zs;
	++tricount;

	// bottom
	tri_x[tricount][0] = xs;
	tri_y[tricount][0] = ys;
	tri_z[tricount][0] = ze;

	tri_x[tricount][1] = xe;
	tri_y[tricount][1] = ys;
	tri_z[tricount][1] = ze;

	tri_x[tricount][2] = xe;
	tri_y[tricount][2] = ye;
	tri_z[tricount][2] = ze;
	++tricount;

	tri_x[tricount][0] = xs;
	tri_y[tricount][0] = ys;
	tri_z[tricount][0] = ze;

	tri_x[tricount][1] = xe;
	tri_y[tricount][1] = ye;
	tri_z[tricount][1] = ze;

	tri_x[tricount][2] = xs;
	tri_y[tricount][2] = ye;
	tri_z[tricount][2] = ze;
	++tricount;

	// front
	tri_x[tricount][0] = xs;
	tri_y[tricount][0] = ys;
	tri_z[tricount][0] = zs;

	tri_x[tricount][1] = xs;
	tri_y[tricount][1] = ys;
	tri_z[tricount][1] = ze;

	tri_x[tricount][2] = xs;
	tri_y[tricount][2] = ye;
	tri_z[tricount][2] = ze;
	++tricount;

	tri_x[tricount][0] = xs;
	tri_y[tricount][0] = ye;
	tri_z[tricount][0] = ze;

	tri_x[tricount][1] = xs;
	tri_y[tricount][1] = ye;
	tri_z[tricount][1] = zs;

	tri_x[tricount][2] = xs;
	tri_y[tricount][2] = ys;
	tri_z[tricount][2] = zs;
	++tricount;

	// top
	tri_x[tricount][0] = xs;
	tri_y[tricount][0] = ys;
	tri_z[tricount][0] = zs;

	tri_x[tricount][1] = xe;
	tri_y[tricount][1] = ys;
	tri_z[tricount][1] = ze;

	tri_x[tricount][2] = xe;
	tri_y[tricount][2] = ye;
	tri_z[tricount][2] = ze;
	++tricount;

	tri_x[tricount][0] = xe;
	tri_y[tricount][0] = ye;
	tri_z[tricount][0] = ze;

	tri_x[tricount][1] = xs;
	tri_y[tricount][1] = ye;
	tri_z[tricount][1] = zs;

	tri_x[tricount][2] = xs;
	tri_y[tricount][2] = ys;
	tri_z[tricount][2] = zs;
	++tricount;
	}
	 
	tend[entity_count]=tricount;
    
}