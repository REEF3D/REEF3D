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

#include"geo_primitive.h"
#include"increment.h"
#include<cmath>

// Geometry core: triangulated primitives, taken from the 6DOF floating body objects
// (box, cylinder_x/y/z, wedge_sym, wedge, hexahedron) and the NHFLOW solid forcing
// geometry (sphere, wedge_x/y/z). The triangles are appended at tricount.

int geo_primitive::segments(double r, double ds, int nmin)
{
    const double U = 2.0*PI*r;

    return MAX(int(U/ds),nmin);
}

void geo_primitive::box(double **tri_x, double **tri_y, double **tri_z, int &tricount,
                       double xs, double xe, double ys, double ye, double zs, double ze)
{
    
	
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
}

void geo_primitive::cylinder_x(double **tri_x, double **tri_y, double **tri_z, int &tricount,
                              double ym, double zm, double x1, double x2, double r, int snum)
{
    int n;
    double ds,phi;
    
    ds = (2.0*PI)/double(snum);
    
    phi=0.0;
    
	
	for(n=0;n<snum;++n)
	{
	//bottom circle	
	tri_x[tricount][0] = x1;
	tri_y[tricount][0] = ym;
	tri_z[tricount][0] = zm;
	
	tri_x[tricount][1] = x1;
	tri_y[tricount][1] = ym + r*sin(phi);
	tri_z[tricount][1] = zm + r*cos(phi);
	
	tri_x[tricount][2] = x1;
	tri_y[tricount][2] = ym + r*sin(phi+ds);
	tri_z[tricount][2] = zm + r*cos(phi+ds);
	++tricount;
		
	//top circle
	tri_x[tricount][0] = x2;
	tri_y[tricount][0] = ym;
	tri_z[tricount][0] = zm;
	
	tri_x[tricount][1] = x2;
	tri_y[tricount][1] = ym + r*sin(phi);
	tri_z[tricount][1] = zm + r*cos(phi);
	
	tri_x[tricount][2] = x2;
	tri_y[tricount][2] = ym + r*sin(phi+ds);
	tri_z[tricount][2] = zm + r*cos(phi+ds);
	++tricount;
	
	//side		
	// 1st triangle
	tri_x[tricount][0] = x1;
	tri_y[tricount][0] = ym + r*sin(phi);
	tri_z[tricount][0] = zm + r*cos(phi);
	
	tri_x[tricount][1] = x2;
	tri_y[tricount][1] = ym + r*sin(phi+ds);
	tri_z[tricount][1] = zm + r*cos(phi+ds);
	
	tri_x[tricount][2] = x1;
	tri_y[tricount][2] = ym + r*sin(phi+ds);
	tri_z[tricount][2] = zm + r*cos(phi+ds);

	++tricount;
	
	// 2nd triangle
	tri_x[tricount][0] = x1;
	tri_y[tricount][0] = ym + r*sin(phi);
	tri_z[tricount][0] = zm + r*cos(phi);
	
	tri_x[tricount][1] = x2;
	tri_y[tricount][1] = ym + r*sin(phi+ds);
	tri_z[tricount][1] = zm + r*cos(phi+ds);
	
	tri_x[tricount][2] = x2;
	tri_y[tricount][2] = ym + r*sin(phi);
	tri_z[tricount][2] = zm + r*cos(phi);
	
	++tricount;
		
	phi+=ds;
	}
}

void geo_primitive::cylinder_y(double **tri_x, double **tri_y, double **tri_z, int &tricount,
                              double xm, double zm, double y1, double y2, double r, int snum)
{
    int n;
    double ds,phi;
    
    ds = (2.0*PI)/double(snum);
    
    phi=0.0;
    

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
}

void geo_primitive::cylinder_z(double **tri_x, double **tri_y, double **tri_z, int &tricount,
                              double xm, double ym, double z1, double z2, double r, int snum)
{
    int n;
    double ds,phi;
    
    ds = (2.0*PI)/double(snum);
    
    phi=0.0;
    
	
	for(n=0;n<snum;++n)
	{
	//bottom circle	
	tri_x[tricount][0] = xm;
	tri_y[tricount][0] = ym;
	tri_z[tricount][0] = z1;
	
	tri_x[tricount][1] = xm + r*cos(phi);
	tri_y[tricount][1] = ym + r*sin(phi);
	tri_z[tricount][1] = z1;
	
	tri_x[tricount][2] = xm + r*cos(phi+ds);
	tri_y[tricount][2] = ym + r*sin(phi+ds);
	tri_z[tricount][2] = z1;
	++tricount;
		
	//top circle
	tri_x[tricount][0] = xm;
	tri_y[tricount][0] = ym;
	tri_z[tricount][0] = z2;
	
	tri_x[tricount][1] = xm + r*cos(phi);
	tri_y[tricount][1] = ym + r*sin(phi);
	tri_z[tricount][1] = z2;
	
	tri_x[tricount][2] = xm + r*cos(phi+ds);
	tri_y[tricount][2] = ym + r*sin(phi+ds);
	tri_z[tricount][2] = z2;
	++tricount;
	
	//side		
	// 1st triangle
	tri_x[tricount][0] = xm + r*cos(phi);
	tri_y[tricount][0] = ym + r*sin(phi);
	tri_z[tricount][0] = z1;
	
	tri_x[tricount][1] = xm + r*cos(phi+ds);
	tri_y[tricount][1] = ym + r*sin(phi+ds);
	tri_z[tricount][1] = z2;
	
	tri_x[tricount][2] = xm + r*cos(phi+ds);
	tri_y[tricount][2] = ym + r*sin(phi+ds);
	tri_z[tricount][2] = z1;

	++tricount;
	
	// 2nd triangle
	tri_x[tricount][0] = xm + r*cos(phi);
	tri_y[tricount][0] = ym + r*sin(phi);
	tri_z[tricount][0] = z1;
	
	tri_x[tricount][1] = xm + r*cos(phi+ds);
	tri_y[tricount][1] = ym + r*sin(phi+ds);
	tri_z[tricount][1] = z2;
	
	tri_x[tricount][2] = xm + r*cos(phi);
	tri_y[tricount][2] = ym + r*sin(phi);
	tri_z[tricount][2] = z2;
	
	++tricount;
		
	phi+=ds;
	}
}

void geo_primitive::wedge_sym(double **tri_x, double **tri_y, double **tri_z, int &tricount,
                             double xs, double xe, double ys, double ye, double zs, double ze)
{
    double xm;
    
    xm = xs + 0.5*(xe-xs);
    
	
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
}

void geo_primitive::wedge(double **tri_x, double **tri_y, double **tri_z, int &tricount, const double *v)
{
    const double x1=v[0], y1=v[1], z1=v[2], x2=v[3], y2=v[4], z2=v[5], x3=v[6], y3=v[7], z3=v[8];
    const double x4=v[9], y4=v[10], z4=v[11], x5=v[12], y5=v[13], z5=v[14], x6=v[15], y6=v[16], z6=v[17];
    
	
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
}

void geo_primitive::hexahedron(double **tri_x, double **tri_y, double **tri_z, int &tricount, const double *v)
{
    const double x1=v[0], y1=v[1], z1=v[2], x2=v[3], y2=v[4], z2=v[5], x3=v[6], y3=v[7], z3=v[8];
    const double x4=v[9], y4=v[10], z4=v[11], x5=v[12], y5=v[13], z5=v[14], x6=v[15], y6=v[16], z6=v[17];
    const double x7=v[18], y7=v[19], z7=v[20], x8=v[21], y8=v[22], z8=v[23];
    
	
	// Face 3
	// Tri 1
	tri_x[tricount][0] = x1;
	tri_y[tricount][0] = y1;
	tri_z[tricount][0] = z1;
	
	tri_x[tricount][1] = x2;
	tri_y[tricount][1] = y2;
	tri_z[tricount][1] = z2;
	
	tri_x[tricount][2] = x6;
	tri_y[tricount][2] = y6;
	tri_z[tricount][2] = z6;
	++tricount;

	// Tri 2
	tri_x[tricount][0] = x5;
	tri_y[tricount][0] = y5;
	tri_z[tricount][0] = z5;

	tri_x[tricount][1] = x1;
	tri_y[tricount][1] = y1;
	tri_z[tricount][1] = z1;
	
	tri_x[tricount][2] = x6;
	tri_y[tricount][2] = y6;
	tri_z[tricount][2] = z6;
	++tricount;
    
	// Face 4
    // Tri 3
	tri_x[tricount][0] = x7;
	tri_y[tricount][0] = y7;
	tri_z[tricount][0] = z7;
	
	tri_x[tricount][1] = x6;
	tri_y[tricount][1] = y6;
	tri_z[tricount][1] = z6;
	
	tri_x[tricount][2] = x2;
	tri_y[tricount][2] = y2;
	tri_z[tricount][2] = z2;
	++tricount;
	
	// Tri 4
	tri_x[tricount][0] = x3;
	tri_y[tricount][0] = y3;
	tri_z[tricount][0] = z3;
	
	tri_x[tricount][1] = x7;
	tri_y[tricount][1] = y7;
	tri_z[tricount][1] = z7;
	
	tri_x[tricount][2] = x2;
	tri_y[tricount][2] = y2;
	tri_z[tricount][2] = z2;
	++tricount;
	
	// Face 1
	// Tri 5
	tri_x[tricount][0] = x8;
	tri_y[tricount][0] = y8;
	tri_z[tricount][0] = z8;
	
	tri_x[tricount][1] = x4;
	tri_y[tricount][1] = y4;
	tri_z[tricount][1] = z4;
	
	tri_x[tricount][2] = x1;
	tri_y[tricount][2] = y1;
	tri_z[tricount][2] = z1;
	++tricount;
	
	// Tri 6
	tri_x[tricount][0] = x5;
	tri_y[tricount][0] = y5;
	tri_z[tricount][0] = z5;
	
	tri_x[tricount][1] = x8;
	tri_y[tricount][1] = y8;
	tri_z[tricount][1] = z8;
	
	tri_x[tricount][2] = x1;
	tri_y[tricount][2] = y1;
	tri_z[tricount][2] = z1;
	++tricount;
	
	// Face 2
	// Tri 7
	tri_x[tricount][0] = x7;
	tri_y[tricount][0] = y7;
	tri_z[tricount][0] = z7;

	tri_x[tricount][1] = x3;
	tri_y[tricount][1] = y3;
	tri_z[tricount][1] = z3;
	
	tri_x[tricount][2] = x4;
	tri_y[tricount][2] = y4;
	tri_z[tricount][2] = z4;
	++tricount;
	
	// Tri 8
	tri_x[tricount][0] = x8;
	tri_y[tricount][0] = y8;
	tri_z[tricount][0] = z8;
	
	tri_x[tricount][1] = x7;
	tri_y[tricount][1] = y7;
	tri_z[tricount][1] = z7;
	
	tri_x[tricount][2] = x4;
	tri_y[tricount][2] = y4;
	tri_z[tricount][2] = z4;
	++tricount;
	
	// Face 5
	// Tri 9
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
	
	// Tri 10
	tri_x[tricount][0] = x4;
	tri_y[tricount][0] = y4;
	tri_z[tricount][0] = z4;
	
	tri_x[tricount][1] = x3;
	tri_y[tricount][1] = y3;
	tri_z[tricount][1] = z3;
	
	tri_x[tricount][2] = x1;
	tri_y[tricount][2] = y1;
	tri_z[tricount][2] = z1;
	++tricount;
	
	// Face 6
	// Tri 11
	tri_x[tricount][0] = x5;
	tri_y[tricount][0] = y5;
	tri_z[tricount][0] = z5;
	
	tri_x[tricount][1] = x6;
	tri_y[tricount][1] = y6;
	tri_z[tricount][1] = z6;
	
	tri_x[tricount][2] = x7;
	tri_y[tricount][2] = y7;
	tri_z[tricount][2] = z7;
	++tricount;
	
	// Tri 12
	tri_x[tricount][0] = x5;
	tri_y[tricount][0] = y5;
	tri_z[tricount][0] = z5;
	
	tri_x[tricount][1] = x7;
	tri_y[tricount][1] = y7;
	tri_z[tricount][1] = z7;
	
	tri_x[tricount][2] = x8;
	tri_y[tricount][2] = y8;
	tri_z[tricount][2] = z8;
	++tricount;
}

void geo_primitive::sphere(double **tri_x, double **tri_y, double **tri_z, int &tricount,
                          double xm, double ym, double zm, double r, int snum)
{
    int n,q;
    double ds,dt,phi,theta;
    
    ds = (2.0*PI)/double(snum);
    dt = ds;
    
    phi=-0.5*PI;
    theta=-0.5*PI;
    

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
}

void geo_primitive::wedge_x(double **tri_x, double **tri_y, double **tri_z, int &tricount,
                           double xs, double xe, double ys, double ye, double zs, double ze)
{
    
	
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
}

void geo_primitive::wedge_y(double **tri_x, double **tri_y, double **tri_z, int &tricount,
                           double xs, double xe, double ys, double ye, double zs, double ze)
{
    
	
	if(zs<ze)
	{
	// sides
	tri_x[tricount][0] = xs;
	tri_y[tricount][0] = ys;
	tri_z[tricount][0] = zs;

	tri_x[tricount][1] = xs;
	tri_y[tricount][1] = ye;
	tri_z[tricount][1] = zs;

	tri_x[tricount][2] = xs;
	tri_y[tricount][2] = ye;
	tri_z[tricount][2] = ze;
	++tricount;

	tri_x[tricount][0] = xe;
	tri_y[tricount][0] = ys;
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

	tri_x[tricount][1] = xs;
	tri_y[tricount][1] = ye;
	tri_z[tricount][1] = zs;

	tri_x[tricount][2] = xe;
	tri_y[tricount][2] = ye;
	tri_z[tricount][2] = zs;
	++tricount;

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

	// front
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

	tri_x[tricount][0] = xs;
	tri_y[tricount][0] = ye;
	tri_z[tricount][0] = zs;

	tri_x[tricount][1] = xs;
	tri_y[tricount][1] = ye;
	tri_z[tricount][1] = ze;

	tri_x[tricount][2] = xe;
	tri_y[tricount][2] = ye;
	tri_z[tricount][2] = ze;
	++tricount;

	// top
	tri_x[tricount][0] = xs;
	tri_y[tricount][0] = ys;
	tri_z[tricount][0] = zs;

	tri_x[tricount][1] = xe;
	tri_y[tricount][1] = ye;
	tri_z[tricount][1] = ze;

	tri_x[tricount][2] = xs;
	tri_y[tricount][2] = ye;
	tri_z[tricount][2] = ze;
	++tricount;

	tri_x[tricount][0] = xs;
	tri_y[tricount][0] = ys;
	tri_z[tricount][0] = zs;

	tri_x[tricount][1] = xe;
	tri_y[tricount][1] = ys;
	tri_z[tricount][1] = zs;

	tri_x[tricount][2] = xe;
	tri_y[tricount][2] = ye;
	tri_z[tricount][2] = ze;
	++tricount;
	}

	if(zs>=ze)
	{
	// sides
	tri_x[tricount][0] = xs;
	tri_y[tricount][0] = ys;
	tri_z[tricount][0] = ze;

	tri_x[tricount][1] = xs;
	tri_y[tricount][1] = ye;
	tri_z[tricount][1] = ze;

	tri_x[tricount][2] = xs;
	tri_y[tricount][2] = ys;
	tri_z[tricount][2] = zs;
	++tricount;

	tri_x[tricount][0] = xe;
	tri_y[tricount][0] = ys;
	tri_z[tricount][0] = ze;

	tri_x[tricount][1] = xe;
	tri_y[tricount][1] = ye;
	tri_z[tricount][1] = ze;

	tri_x[tricount][2] = xe;
	tri_y[tricount][2] = ys;
	tri_z[tricount][2] = zs;
	++tricount;

	// bottom
	tri_x[tricount][0] = xs;
	tri_y[tricount][0] = ys;
	tri_z[tricount][0] = ze;

	tri_x[tricount][1] = xs;
	tri_y[tricount][1] = ye;
	tri_z[tricount][1] = ze;

	tri_x[tricount][2] = xe;
	tri_y[tricount][2] = ye;
	tri_z[tricount][2] = ze;
	++tricount;

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

	// front
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
	tri_y[tricount][0] = ys;
	tri_z[tricount][0] = zs;

	tri_x[tricount][1] = xs;
	tri_y[tricount][1] = ys;
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
	tri_y[tricount][1] = ye;
	tri_z[tricount][1] = ze;

	tri_x[tricount][2] = xs;
	tri_y[tricount][2] = ye;
	tri_z[tricount][2] = ze;
	++tricount;
    
	tri_x[tricount][0] = xs;
	tri_y[tricount][0] = ys;
	tri_z[tricount][0] = zs;

	tri_x[tricount][1] = xe;
	tri_y[tricount][1] = ys;
	tri_z[tricount][1] = zs;

	tri_x[tricount][2] = xe;
	tri_y[tricount][2] = ye;
	tri_z[tricount][2] = ze;
	++tricount;
	}
}

void geo_primitive::wedge_z(double **tri_x, double **tri_y, double **tri_z, int &tricount,
                           double xs, double xe, double ys, double ye, double zs, double ze)
{
    
	
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
}

