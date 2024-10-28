/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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
#include"fdm_nhf.h"
#include"ghostcell.h"

void nhflow_forcing::cylinder_y(lexer *p, ghostcell *pgc, int id)
{
    double U,ds,phi;
	double r1,y1,y2;
    double xm,zm;
	int snum;
	
	xm=p->A583_xc[id];
    zm=p->A583_zc[id];
	
	y1=p->A583_ys[id];
	y2=p->A583_ye[id];
	
    r1=p->A583_r[id];

	U = 2.0 * PI * r1;	
	ds = 0.75*(U*p->DXM);
	snum = int(U/ds);

// Vertices	
	ds = (2.0*PI)/double(snum);
	
	phi=0.0;
	
	tstart[entity_count]=tricount;
	

	for(n=0;n<snum;++n)
	{
	//bottom circle	
	tri_x[tricount][0] = xm;
	tri_y[tricount][0] = y1;
	tri_z[tricount][0] = zm;
	
	tri_x[tricount][1] = xm + r1*sin(phi);
	tri_y[tricount][1] = y1;
	tri_z[tricount][1] = zm + r1*cos(phi);
	
	tri_x[tricount][2] = xm + r1*sin(phi+ds);
	tri_y[tricount][2] = y1;
	tri_z[tricount][2] = zm + r1*cos(phi+ds);
	++tricount;
		
	//top circle
	tri_x[tricount][0] = xm;
	tri_y[tricount][0] = y2;
	tri_z[tricount][0] = zm;
	
	tri_x[tricount][1] = xm + r1*sin(phi);
	tri_y[tricount][1] = y2;
	tri_z[tricount][1] = zm + r1*cos(phi);
	
	tri_x[tricount][2] = xm + r1*sin(phi+ds);
	tri_y[tricount][2] = y2;
	tri_z[tricount][2] = zm + r1*cos(phi+ds);
	++tricount;
	
	//side		
	// 1st triangle
	tri_x[tricount][0] = xm + r1*sin(phi);
	tri_y[tricount][0] = y1;
	tri_z[tricount][0] = zm + r1*cos(phi);
	
	tri_x[tricount][1] = xm + r1*sin(phi+ds);
	tri_y[tricount][1] = y2;
	tri_z[tricount][1] = zm + r1*cos(phi+ds);
	
	tri_x[tricount][2] = xm + r1*sin(phi+ds);
	tri_y[tricount][2] = y1;
	tri_z[tricount][2] = zm + r1*cos(phi+ds);

	++tricount;
	
	// 2nd triangle
	tri_x[tricount][0] = xm + r1*sin(phi);
	tri_y[tricount][0] = y1;
	tri_z[tricount][0] = zm + r1*cos(phi);
	
	tri_x[tricount][1] = xm + r1*sin(phi+ds);
	tri_y[tricount][1] = y2;
	tri_z[tricount][1] = zm + r1*cos(phi+ds);
	
	tri_x[tricount][2] = xm + r1*sin(phi);
	tri_y[tricount][2] = y2;
	tri_z[tricount][2] = zm + r1*cos(phi);
	++tricount;

	phi+=ds;
	}
	
	tend[entity_count]=tricount;
    
    
}