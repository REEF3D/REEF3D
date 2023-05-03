/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

void sixdof_df_object::wedge(lexer *p, fdm *a, ghostcell *pgc, int id)
{

    double x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;
    double x5,y5,z5,x6,y6,z6;
	
	
	x1 = p->X163_x1[id];
    y1 = p->X163_y1[id];
    z1 = p->X163_z1[id];
    
    x2 = p->X163_x2[id];
    y2 = p->X163_y2[id];
    z2 = p->X163_z2[id];
    
    x3 = p->X163_x3[id];
    y3 = p->X163_y3[id];
    z3 = p->X163_z3[id];
    
    x4 = p->X163_x4[id];
    y4 = p->X163_y4[id];
    z4 = p->X163_z4[id];
	
	x5 = p->X163_x5[id];
    y5 = p->X163_y5[id];
    z5 = p->X163_z5[id];
	
	x6 = p->X163_x6[id];
    y6 = p->X163_y6[id];
    z6 = p->X163_z6[id];
  
	
	tstart[entity_count]=tricount;
	
	tri_x[tricount][0] = x3;
	tri_y[tricount][0] = y3;
	tri_z[tricount][0] = z3;
	
	tri_x[tricount][1] = x2;
	tri_y[tricount][1] = y2;
	tri_z[tricount][1] = z2;
	
	tri_x[tricount][2] = x1;
	tri_y[tricount][2] = y1;
	tri_z[tricount][2] = z1;
	++tricount;
	
	tri_x[tricount][0] = x1;
	tri_y[tricount][0] = y1;
	tri_z[tricount][0] = z1;

	tri_x[tricount][1] = x2;
	tri_y[tricount][1] = y2;
	tri_z[tricount][1] = z2;
	
	tri_x[tricount][2] = x5;
	tri_y[tricount][2] = y5;
	tri_z[tricount][2] = z5;
	++tricount;
    
	
	tri_x[tricount][0] = x5;
	tri_y[tricount][0] = y5;
	tri_z[tricount][0] = z5;
	
	tri_x[tricount][1] = x4;
	tri_y[tricount][1] = y4;
	tri_z[tricount][1] = z4;
	
	tri_x[tricount][2] = x1;
	tri_y[tricount][2] = y1;
	tri_z[tricount][2] = z1;
	++tricount;

	
	tri_x[tricount][0] = x2;
	tri_y[tricount][0] = y2;
	tri_z[tricount][0] = z2;
	
	tri_x[tricount][1] = x3;
	tri_y[tricount][1] = y3;
	tri_z[tricount][1] = z3;
	
	tri_x[tricount][2] = x5;
	tri_y[tricount][2] = y5;
	tri_z[tricount][2] = z5;
	++tricount;
	
	
	tri_x[tricount][0] = x6;
	tri_y[tricount][0] = y6;
	tri_z[tricount][0] = z6;
	
	tri_x[tricount][1] = x5;
	tri_y[tricount][1] = y5;
	tri_z[tricount][1] = z5;
	
	tri_x[tricount][2] = x3;
	tri_y[tricount][2] = y3;
	tri_z[tricount][2] = z3;
	++tricount;
	
	
	tri_x[tricount][0] = x6;
	tri_y[tricount][0] = y6;
	tri_z[tricount][0] = z6;
	
	tri_x[tricount][1] = x3;
	tri_y[tricount][1] = y3;
	tri_z[tricount][1] = z3;
	
	tri_x[tricount][2] = x1;
	tri_y[tricount][2] = y1;
	tri_z[tricount][2] = z1;
	++tricount;
	
	
	tri_x[tricount][0] = x4;
	tri_y[tricount][0] = y4;
	tri_z[tricount][0] = z4;
	
	tri_x[tricount][1] = x6;
	tri_y[tricount][1] = y6;
	tri_z[tricount][1] = z6;
	
	tri_x[tricount][2] = x1;
	tri_y[tricount][2] = y1;
	tri_z[tricount][2] = z1;
	++tricount;
	
	
	tri_x[tricount][0] = x4;
	tri_y[tricount][0] = y4;
	tri_z[tricount][0] = z4;
	
	tri_x[tricount][1] = x5;
	tri_y[tricount][1] = y5;
	tri_z[tricount][1] = z5;
	
	tri_x[tricount][2] = x6;
	tri_y[tricount][2] = y6;
	tri_z[tricount][2] = z6;
	++tricount;
	
	
	tend[entity_count]=tricount;
    
}

