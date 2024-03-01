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

#include"6DOF_obj.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void sixdof_obj::wedge_sym(lexer *p, fdm *a, ghostcell *pgc, int id)
{
	double xm;
		
	
	xs = p->X153_xs;
    xe = p->X153_xe;
	
    ys = p->X153_ys;
    ye = p->X153_ye;

    zs = p->X153_zs;
    ze = p->X153_ze;  

	xm = xs + 0.5*(xe-xs);
	

// Vert

// Face 3
	// Tri 1
	tstart[entity_count]=tricount;
	
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
	
	// Tri 2
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

// Sides
	// Tri 3
	tri_x[tricount][0] = xs;
	tri_y[tricount][0] = ys;
	tri_z[tricount][0] = zs;
	
	tri_x[tricount][1] = xe;
	tri_y[tricount][1] = ys;
	tri_z[tricount][1] = zs;
	
	tri_x[tricount][2] = xm;
	tri_y[tricount][2] = ys;
	tri_z[tricount][2] = ze;
	++tricount;
	
	// Tri 4
	tri_x[tricount][0] = xs;
	tri_y[tricount][0] = ye;
	tri_z[tricount][0] = zs;
	
	tri_x[tricount][1] = xe;
	tri_y[tricount][1] = ye;
	tri_z[tricount][1] = zs;
	
	tri_x[tricount][2] = xm;
	tri_y[tricount][2] = ye;
	tri_z[tricount][2] = ze;
	++tricount;

// Front	
	// Tri 5
	tri_x[tricount][0] = xs;
	tri_y[tricount][0] = ys;
	tri_z[tricount][0] = zs;
	
	tri_x[tricount][1] = xm;
	tri_y[tricount][1] = ys;
	tri_z[tricount][1] = ze;
	
	tri_x[tricount][2] = xs;
	tri_y[tricount][2] = ye;
	tri_z[tricount][2] = zs;
	++tricount;
	
	// Tri 6
	tri_x[tricount][0] = xs;
	tri_y[tricount][0] = ye;
	tri_z[tricount][0] = zs;
	
	tri_x[tricount][1] = xm;
	tri_y[tricount][1] = ye;
	tri_z[tricount][1] = ze;
	
	tri_x[tricount][2] = xm;
	tri_y[tricount][2] = ys;
	tri_z[tricount][2] = ze;
	++tricount;

// Back	
	// Tri 7
	tri_x[tricount][0] = xm;
	tri_y[tricount][0] = ys;
	tri_z[tricount][0] = ze;
	
	tri_x[tricount][1] = xe;
	tri_y[tricount][1] = ys;
	tri_z[tricount][1] = zs;
	
	tri_x[tricount][2] = xe;
	tri_y[tricount][2] = ye;
	tri_z[tricount][2] = zs;
	++tricount;
	
	// Tri 8	
	tri_x[tricount][0] = xm;
	tri_y[tricount][0] = ys;
	tri_z[tricount][0] = ze;
	
	tri_x[tricount][1] = xm;
	tri_y[tricount][1] = ye;
	tri_z[tricount][1] = ze;
	
	tri_x[tricount][2] = xe;
	tri_y[tricount][2] = ye;
	tri_z[tricount][2] = zs;
	++tricount;

	
	tend[entity_count]=tricount;
}

