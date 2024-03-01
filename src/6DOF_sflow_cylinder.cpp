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
--------------------------------------------------------------------*/

#include"6DOF_sflow.h"
#include"lexer.h"
#include"ghostcell.h"

void sixdof_sflow::cylinder(lexer *p, ghostcell *pgc)
{
	double U,ds,angle;
	double xm,ym,zm,z1,z2,r;
	int snum, trisum;
	
	xm = p->X133_xc;
	ym = p->X133_yc;
    zm = p->X133_zc;
	z1 = zm - 0.5*p->X133_h;
	z2 = zm + 0.5*p->X133_h;
    r = p->X133_rad;

    // Prepare fields
    U = 2.0 * PI * r;
	ds = p->DXM;
	snum = int(U/ds);
    trisum=5*(snum+1);
    p->Darray(tri_x,trisum,3);
	p->Darray(tri_y,trisum,3);
	p->Darray(tri_z,trisum,3);
    p->Darray(tri_x0,trisum,3);
	p->Darray(tri_y0,trisum,3);
	p->Darray(tri_z0,trisum,3);    	

    // Vertices	
	ds = (2.0*PI)/double(snum);
	angle=0.0;
    tricount = 0;

	for(int n=0; n<snum; ++n)
	{
        //bottom circle	
        tri_x[tricount][0] = xm;
        tri_y[tricount][0] = ym;
        tri_z[tricount][0] = z1;
        
        tri_x[tricount][1] = xm + r*cos(angle);
        tri_y[tricount][1] = ym + r*sin(angle);
        tri_z[tricount][1] = z1;
        
        tri_x[tricount][2] = xm + r*cos(angle+ds);
        tri_y[tricount][2] = ym + r*sin(angle+ds);
        tri_z[tricount][2] = z1;
        ++tricount;
            
        //top circle
        tri_x[tricount][0] = xm;
        tri_y[tricount][0] = ym;
        tri_z[tricount][0] = z2;
        
        tri_x[tricount][1] = xm + r*cos(angle);
        tri_y[tricount][1] = ym + r*sin(angle);
        tri_z[tricount][1] = z2;
        
        tri_x[tricount][2] = xm + r*cos(angle+ds);
        tri_y[tricount][2] = ym + r*sin(angle+ds);
        tri_z[tricount][2] = z2;
        ++tricount;
        
        //side		
        // 1st triangle
        tri_x[tricount][0] = xm + r*cos(angle);
        tri_y[tricount][0] = ym + r*sin(angle);
        tri_z[tricount][0] = z1;
        
        tri_x[tricount][1] = xm + r*cos(angle+ds);
        tri_y[tricount][1] = ym + r*sin(angle+ds);
        tri_z[tricount][1] = z2;
        
        tri_x[tricount][2] = xm + r*cos(angle+ds);
        tri_y[tricount][2] = ym + r*sin(angle+ds);
        tri_z[tricount][2] = z1;
        ++tricount;
        
        // 2nd triangle
        tri_x[tricount][0] = xm + r*cos(angle);
        tri_y[tricount][0] = ym + r*sin(angle);
        tri_z[tricount][0] = z1;
        
        tri_x[tricount][1] = xm + r*cos(angle+ds);
        tri_y[tricount][1] = ym + r*sin(angle+ds);
        tri_z[tricount][1] = z2;
        
        tri_x[tricount][2] = xm + r*cos(angle);
        tri_y[tricount][2] = ym + r*sin(angle);
        tri_z[tricount][2] = z2;
        ++tricount;
            
        angle+=ds;
	}
}

