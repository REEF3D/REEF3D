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
#include"fdm_nhf.h"
#include"ghostcell.h"

void nhflow_forcing::cylinder_z(lexer *p, ghostcell *pgc, int id)
{
    double U,ds,phi;
	double r1,z1,z2;
    double xm,ym;
	int snum;
	
    
	xm=p->A584_xc[id];
    ym=p->A584_yc[id];
	
	z1=p->A584_zs[id];
	z2=p->A584_ze[id];
	
    r1=p->A584_r[id];
    
    //cout<<p->mpirank<<" FORCING_CYLINDER  xm: "<<xm<<" ym: "<<ym<<" r1: "<<r1<<endl;


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
	tri_y[tricount][0] = ym;
	tri_z[tricount][0] = z1;
	
	tri_x[tricount][1] = xm + r1*cos(phi);
	tri_y[tricount][1] = ym + r1*sin(phi);
	tri_z[tricount][1] = z1;
	
	tri_x[tricount][2] = xm + r1*cos(phi+ds);
	tri_y[tricount][2] = ym + r1*sin(phi+ds);
	tri_z[tricount][2] = z1;
	++tricount;
		
	//top circle
	tri_x[tricount][0] = xm;
	tri_y[tricount][0] = ym;
	tri_z[tricount][0] = z2;
	
	tri_x[tricount][1] = xm + r1*cos(phi);
	tri_y[tricount][1] = ym + r1*sin(phi);
	tri_z[tricount][1] = z2;
	
	tri_x[tricount][2] = xm + r1*cos(phi+ds);
	tri_y[tricount][2] = ym + r1*sin(phi+ds);
	tri_z[tricount][2] = z2;
	++tricount;
	
	//side		
	// 1st triangle
	tri_x[tricount][0] = xm + r1*cos(phi);
	tri_y[tricount][0] = ym + r1*sin(phi);
	tri_z[tricount][0] = z1;
	
	tri_x[tricount][1] = xm + r1*cos(phi+ds);
	tri_y[tricount][1] = ym + r1*sin(phi+ds);
	tri_z[tricount][1] = z2;
	
	tri_x[tricount][2] = xm + r1*cos(phi+ds);
	tri_y[tricount][2] = ym + r1*sin(phi+ds);
	tri_z[tricount][2] = z1;

	++tricount;
	
	// 2nd triangle
	tri_x[tricount][0] = xm + r1*cos(phi);
	tri_y[tricount][0] = ym + r1*sin(phi);
	tri_z[tricount][0] = z1;
	
	tri_x[tricount][1] = xm + r1*cos(phi+ds);
	tri_y[tricount][1] = ym + r1*sin(phi+ds);
	tri_z[tricount][1] = z2;
	
	tri_x[tricount][2] = xm + r1*cos(phi);
	tri_y[tricount][2] = ym + r1*sin(phi);
	tri_z[tricount][2] = z2;
	
	++tricount;
		
	phi+=ds;
	}
		
	
	tend[entity_count]=tricount;
    
}