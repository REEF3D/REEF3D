/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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

#include"6DOF_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void sixdof_f::cylinder_y(lexer *p, fdm *a, ghostcell *pgc, int id)
{
	double U,ds,phi;
	double xm,ym,zm,y1,y2,r;
	int snum;
	int vertice_mem, center1_num,center2_num;
	double Ax,Ay,Az;
	double Bx,By,Bz;
	double Nx,Ny,Nz,norm;
	
	

	xm=p->X132_xc;
	ym=p->X132_yc;
    zm=p->X132_zc;
	
	y1=ym-0.5*p->X132_h;
	y2=ym+0.5*p->X132_h;
	
    r=p->X132_rad;
	
	U = 2.0 * PI * r;
	
	ds = 0.25*(U*p->dx);
	
	snum = int(U/ds);
	
	

// Vertices	
	ds = (2.0*PI)/double(snum);
	
	phi=0.0;

	tstart[entity_count]=tricount;
	
	//cout<<p->mpirank<<" DISC: "<<tricount<<" . "<<xm<<" "<<ym<<" "<<zm<<" "<<y1<<" "<<y2<<" "<<r<<" "<<endl;
	
	for(n=0;n<snum;++n)
	{
	//bottom circle	
	tri_x[tricount][0] = xm;
	tri_y[tricount][0] = y1;
	tri_z[tricount][0] = zm;
	
	tri_x[tricount][1] = xm + r*sin(phi);
	tri_y[tricount][1] = y1;
	tri_z[tricount][1] = zm + r*cos(phi);
	
	tri_x[tricount][2] = xm + r*sin(phi+ds);
	tri_y[tricount][2] = y1;
	tri_z[tricount][2] = zm + r*cos(phi+ds);
	
	
	
	++tricount;
		
		
		
	//top circle
	tri_x[tricount][0] = xm;
	tri_y[tricount][0] = y2;
	tri_z[tricount][0] = zm;
	
	tri_x[tricount][1] = xm + r*sin(phi);
	tri_y[tricount][1] = y2;
	tri_z[tricount][1] = zm + r*cos(phi);
	
	tri_x[tricount][2] = xm + r*sin(phi+ds);
	tri_y[tricount][2] = y2;
	tri_z[tricount][2] = zm + r*cos(phi+ds);
	++tricount;
	
	//side		
	// 1st triangle
	tri_x[tricount][0] = xm + r*sin(phi);
	tri_y[tricount][0] = y1;
	tri_z[tricount][0] = zm + r*cos(phi);
	
	tri_x[tricount][1] = xm + r*sin(phi+ds);
	tri_y[tricount][1] = y2;
	tri_z[tricount][1] = zm + r*cos(phi+ds);
	
	tri_x[tricount][2] = xm + r*sin(phi+ds);
	tri_y[tricount][2] = y1;
	tri_z[tricount][2] = zm + r*cos(phi+ds);

	++tricount;
	
	// 2nd triangle
	tri_x[tricount][0] = xm + r*sin(phi);
	tri_y[tricount][0] = y1;
	tri_z[tricount][0] = zm + r*cos(phi);
	
	tri_x[tricount][1] = xm + r*sin(phi+ds);
	tri_y[tricount][1] = y2;
	tri_z[tricount][1] = zm + r*cos(phi+ds);
	
	tri_x[tricount][2] = xm + r*sin(phi);
	tri_y[tricount][2] = y2;
	tri_z[tricount][2] = zm + r*cos(phi);
	++tricount;
	
	
		
	phi+=ds;
	}
		
	tend[entity_count]=tricount;
	
}

