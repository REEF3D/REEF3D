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
#include"ghostcell.h"

void nhflow_forcing::sphere(lexer *p, ghostcell *pgc, int id)
{
    
    double U,ds,dt,phi,theta;
	
	int snum;
	int vertice_mem, center1_num,center2_num;
    double xm,ym,zm,r;
    int q;


	xm=p->A585_xm[id];
    ym=p->A585_ym[id];
    zm=p->A585_zm[id];
    r=p->A585_r[id];


	U = 2.0 * PI * r;
	ds = p->S19*(U*p->DXM);
	snum = int(U/ds);


// Vertices
	ds = (2.0*PI)/double(snum);

    dt = ds;

	phi=-0.5*PI;
	theta=-0.5*PI;

	tstart[entity_count]=tricount;

        // bottom start /triangles
        for(q=0;q<snum;++q)
        {
        tri_x[tricount][0] = xm;
        tri_y[tricount][0] = ym;
        tri_z[tricount][0] = zm-r;

        tri_x[tricount][1] = xm + r*cos(theta+dt)*cos(phi);
        tri_y[tricount][1] = ym + r*cos(theta+dt)*sin(phi);
        tri_z[tricount][1] = zm + r*sin(theta+dt);

        tri_x[tricount][2] = xm + r*cos(theta+dt)*cos(phi+ds);
        tri_y[tricount][2] = ym + r*cos(theta+dt)*sin(phi+ds);
        tri_z[tricount][2] = zm + r*sin(theta+dt);

        ++tricount;

        phi+=ds;
        }

    theta+=dt;

    // middle section / hexahedrons
	for(n=1;n<snum/2-1;++n)
    {
        phi=-0.5*PI;
        for(q=0;q<snum;++q)
        {
        //side
        // 1st triangle
        tri_x[tricount][0] = xm + r*cos(theta)*cos(phi);
        tri_y[tricount][0] = ym + r*cos(theta)*sin(phi);
        tri_z[tricount][0] = zm + r*sin(theta);

        tri_x[tricount][1] = xm + r*cos(theta+dt)*cos(phi);
        tri_y[tricount][1] = ym + r*cos(theta+dt)*sin(phi);
        tri_z[tricount][1] = zm + r*sin(theta+dt);

        tri_x[tricount][2] = xm + r*cos(theta+dt)*cos(phi+ds);
        tri_y[tricount][2] = ym + r*cos(theta+dt)*sin(phi+ds);
        tri_z[tricount][2] = zm + r*sin(theta+dt);

        ++tricount;

        // 2nd triangle
        tri_x[tricount][0] = xm + r*cos(theta)*cos(phi);
        tri_y[tricount][0] = ym + r*cos(theta)*sin(phi);
        tri_z[tricount][0] = zm + r*sin(theta);

        tri_x[tricount][1] = xm + r*cos(theta+dt)*cos(phi+ds);
        tri_y[tricount][1] = ym + r*cos(theta+dt)*sin(phi+ds);
        tri_z[tricount][1] = zm + r*sin(theta+dt);

        tri_x[tricount][2] = xm + r*cos(theta)*cos(phi+ds);
        tri_y[tricount][2] = ym + r*cos(theta)*sin(phi+ds);
        tri_z[tricount][2] = zm + r*sin(theta);

        ++tricount;

        phi+=ds;
        }
    theta+=dt;
	}

    // top start /triangles

        phi=-0.5*PI;
        theta=0.5*PI-dt;
        for(q=0;q<snum;++q)
        {
        tri_x[tricount][0] = xm;
        tri_y[tricount][0] = ym;
        tri_z[tricount][0] = zm+r;

        tri_x[tricount][1] = xm + r*cos(theta)*cos(phi);
        tri_y[tricount][1] = ym + r*cos(theta)*sin(phi);
        tri_z[tricount][1] = zm + r*sin(theta);

        tri_x[tricount][2] = xm + r*cos(theta)*cos(phi+ds);
        tri_y[tricount][2] = ym + r*cos(theta)*sin(phi+ds);
        tri_z[tricount][2] = zm + r*sin(theta);

        ++tricount;

        phi+=ds;
        }

    // end point


	tend[entity_count]=tricount;
}