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

#include"print_porous.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"


void print_porous::box(lexer *p, fdm *a, ghostcell *pgc, int rank)
{
	double xs,ys,zs,xe,ye,ze;
	int vertice_start=vertice_num;
	

	xs = p->B270_xs[rank];
    xe = p->B270_xe[rank];

    ys = p->B270_ys[rank];
    ye = p->B270_ye[rank];

    zs = p->B270_zs[rank];
    ze = p->B270_ze[rank];    
	

	// Vertices	
	vertice[vertice_num][0] = xs;
	vertice[vertice_num][1] = ys;
	vertice[vertice_num][2] = zs;
	++vertice_num;
	
	vertice[vertice_num][0] = xe;
	vertice[vertice_num][1] = ys;
	vertice[vertice_num][2] = zs;
	++vertice_num;
	
	vertice[vertice_num][0] = xe;
	vertice[vertice_num][1] = ye;
	vertice[vertice_num][2] = zs;
	++vertice_num;
	
	vertice[vertice_num][0] = xs;
	vertice[vertice_num][1] = ye;
	vertice[vertice_num][2] = zs;
	++vertice_num;
	
	vertice[vertice_num][0] = xs;
	vertice[vertice_num][1] = ys;
	vertice[vertice_num][2] = ze;
	++vertice_num;
	
	vertice[vertice_num][0] = xe;
	vertice[vertice_num][1] = ys;
	vertice[vertice_num][2] = ze;
	++vertice_num;
	
	vertice[vertice_num][0] = xe;
	vertice[vertice_num][1] = ye;
	vertice[vertice_num][2] = ze;
	++vertice_num;
	
	vertice[vertice_num][0] = xs;
	vertice[vertice_num][1] = ye;
	vertice[vertice_num][2] = ze;
	++vertice_num;
	
	
		
// Polygon

	polygon[polygon_num][0] = 0 + vertice_start;
	polygon[polygon_num][1] = 1 + vertice_start;
	polygon[polygon_num][2] = 2 + vertice_start;
	polygon[polygon_num][3] = 3 + vertice_start;
	numvert[polygon_num] = 4;
	++polygon_num;
	
	polygon[polygon_num][0] = 0 + vertice_start;
	polygon[polygon_num][1] = 1 + vertice_start;
	polygon[polygon_num][2] = 5 + vertice_start;
	polygon[polygon_num][3] = 4 + vertice_start;
	numvert[polygon_num] = 4;
	++polygon_num;
	
	polygon[polygon_num][0] = 1 + vertice_start;
	polygon[polygon_num][1] = 2 + vertice_start;
	polygon[polygon_num][2] = 6 + vertice_start;
	polygon[polygon_num][3] = 5 + vertice_start;
	numvert[polygon_num] = 4;
	++polygon_num;
	
	polygon[polygon_num][0] = 2 + vertice_start;
	polygon[polygon_num][1] = 3 + vertice_start;
	polygon[polygon_num][2] = 7 + vertice_start;
	polygon[polygon_num][3] = 6 + vertice_start;
	numvert[polygon_num] = 4;
	++polygon_num;
	
	polygon[polygon_num][0] = 3 + vertice_start;
	polygon[polygon_num][1] = 0 + vertice_start;
	polygon[polygon_num][2] = 4 + vertice_start;
	polygon[polygon_num][3] = 7 + vertice_start;
	numvert[polygon_num] = 4;
	++polygon_num;
	
	polygon[polygon_num][0] = 4 + vertice_start;
	polygon[polygon_num][1] = 5 + vertice_start;
	polygon[polygon_num][2] = 6 + vertice_start;
	polygon[polygon_num][3] = 7 + vertice_start;
	numvert[polygon_num] = 4;
	++polygon_num;
}

void print_porous::cylinder_z(lexer *p, fdm *a, ghostcell *pgc, int rank)
{
	int vertice_start=vertice_num;
    
    double U,ds,phi;
	double xm,ym,z1,z2;
	double r;
	int snum;
	
	xm=p->B274_xc[rank];
    ym=p->B274_yc[rank];
	
	z1=p->B274_zs[rank];
	z2=p->B274_ze[rank];
	
    r=p->B274_r[rank]; 
	

	U = 2.0 * PI * r;
	
	ds = 0.2*(U*p->DXM);
	
	snum = int(U/ds);
    
cout<<"snum: "<<snum<<" U: "<<U<<" ds: "<<ds<<" r: "<<r<<endl;
// Vertices	
	ds = (2.0*PI)/double(snum);
	
	phi=0.0;
		
    
    for(n=0;n<snum;++n)
	{
	//bottom circle
	vertice[vertice_num][0] = xm;
	vertice[vertice_num][1] = ym;
	vertice[vertice_num][2] = z1;
    ++vertice_num;
	
	vertice[vertice_num][0] = xm + r*cos(phi);
	vertice[vertice_num][1] = ym + r*sin(phi);
	vertice[vertice_num][2] = z1;
    ++vertice_num;
	
	vertice[vertice_num][0] = xm + r*cos(phi+ds);
	vertice[vertice_num][1] = ym + r*sin(phi+ds);
	vertice[vertice_num][2] = z1;
	++vertice_num;
    
		
	//top circle
	vertice[vertice_num][0] = xm;
	vertice[vertice_num][1] = ym;
	vertice[vertice_num][2] = z2;
    ++vertice_num;
	
	vertice[vertice_num][0] = xm + r*cos(phi);
	vertice[vertice_num][1] = ym + r*sin(phi);
	vertice[vertice_num][2] = z2;
    ++vertice_num;
	
	vertice[vertice_num][0] = xm + r*cos(phi+ds);
	vertice[vertice_num][1] = ym + r*sin(phi+ds);
	vertice[vertice_num][2] = z2;
	++vertice_num;
	
	//side		
	// 1st triangle
	vertice[vertice_num][0] = xm + r*cos(phi);
	vertice[vertice_num][1] = ym + r*sin(phi);
	vertice[vertice_num][2] = z1;
    ++vertice_num;
	
	vertice[vertice_num][0] = xm + r*cos(phi+ds);
	vertice[vertice_num][1] = ym + r*sin(phi+ds);
	vertice[vertice_num][2] = z2;
    ++vertice_num;
	
	vertice[vertice_num][0] = xm + r*cos(phi+ds);
	vertice[vertice_num][1] = ym + r*sin(phi+ds);
	vertice[vertice_num][2] = z1;
	++vertice_num;
	
	// 2nd triangle
	vertice[vertice_num][0] = xm + r*cos(phi);
	vertice[vertice_num][1] = ym + r*sin(phi);
	vertice[vertice_num][2] = z1;
    ++vertice_num;
	
	vertice[vertice_num][0] = xm + r*cos(phi+ds);
	vertice[vertice_num][1] = ym + r*sin(phi+ds);
	vertice[vertice_num][2] = z2;
    ++vertice_num;
	
	vertice[vertice_num][0] = xm + r*cos(phi);
	vertice[vertice_num][1] = ym + r*sin(phi);
	vertice[vertice_num][2] = z2;
	++vertice_num;
		
	phi+=ds;
    
    // Polygon
    polygon[polygon_num][0] = 0 + n*12 + vertice_start;
	polygon[polygon_num][1] = 1 + n*12 + vertice_start;
	polygon[polygon_num][2] = 2 + n*12 + vertice_start;
	numvert[polygon_num] = 3;
	++polygon_num;
    
    // Polygon
    polygon[polygon_num][0] = 3 + n*12 + vertice_start;
	polygon[polygon_num][1] = 4 + n*12 + vertice_start;
	polygon[polygon_num][2] = 5 + n*12 + vertice_start;
	numvert[polygon_num] = 3;
	++polygon_num;
    
    // Polygon
    polygon[polygon_num][0] = 6 + n*12 + vertice_start;
	polygon[polygon_num][1] = 7 + n*12 + vertice_start;
	polygon[polygon_num][2] = 8 + n*12 + vertice_start;
	numvert[polygon_num] = 3;
	++polygon_num;
    
    // Polygon
    polygon[polygon_num][0] = 9 + n*12 + vertice_start;
	polygon[polygon_num][1] = 10 + n*12 + vertice_start;
	polygon[polygon_num][2] = 11 + n*12 + vertice_start;
	numvert[polygon_num] = 3;
	++polygon_num;
	}
    
}

void print_porous::wedge_x(lexer *p, fdm *a, ghostcell *pgc,int rank)
{	
	double xs,ys,zs,xe,ye,ze;
	int vertice_start=vertice_num;
	

	xs = p->B281_xs[rank];
    xe = p->B281_xe[rank];

    ys = p->B281_ys[rank];
    ye = p->B281_ye[rank];

    zs = p->B281_zs[rank];
    ze = p->B281_ze[rank];  

	if(zs<ze)
	{
	// Vertices	
		vertice[vertice_num][0] = xs;
		vertice[vertice_num][1] = ys;
		vertice[vertice_num][2] = zs;
		++vertice_num;
		
		vertice[vertice_num][0] = xe;
		vertice[vertice_num][1] = ys;
		vertice[vertice_num][2] = zs;
		++vertice_num;
		
		vertice[vertice_num][0] = xe;
		vertice[vertice_num][1] = ys;
		vertice[vertice_num][2] = ze;
		++vertice_num;
		
		vertice[vertice_num][0] = xs;
		vertice[vertice_num][1] = ye;
		vertice[vertice_num][2] = zs;
		++vertice_num;
		
		vertice[vertice_num][0] = xe;
		vertice[vertice_num][1] = ye;
		vertice[vertice_num][2] = zs;
		++vertice_num;
		
		vertice[vertice_num][0] = xe;
		vertice[vertice_num][1] = ye;
		vertice[vertice_num][2] = ze;
		++vertice_num;
		

	}
	
	if(zs>=ze)
	{
	// Vertices	
		vertice[vertice_num][0] = xs;
		vertice[vertice_num][1] = ys;
		vertice[vertice_num][2] = ze;
		++vertice_num;
		
		vertice[vertice_num][0] = xe;
		vertice[vertice_num][1] = ys;
		vertice[vertice_num][2] = ze;
		++vertice_num;
		
		vertice[vertice_num][0] = xs;
		vertice[vertice_num][1] = ys;
		vertice[vertice_num][2] = zs;
		++vertice_num;
		
		vertice[vertice_num][0] = xs;
		vertice[vertice_num][1] = ye;
		vertice[vertice_num][2] = ze;
		++vertice_num;
		
		vertice[vertice_num][0] = xe;
		vertice[vertice_num][1] = ye;
		vertice[vertice_num][2] = ze;
		++vertice_num;
		
		vertice[vertice_num][0] = xs;
		vertice[vertice_num][1] = ye;
		vertice[vertice_num][2] = zs;
		++vertice_num;
	}
	
		
// Polygon
	if(zs<ze)
	{
	polygon[polygon_num][0] = 0 + vertice_start;
	polygon[polygon_num][1] = 1 + vertice_start;
	polygon[polygon_num][2] = 2 + vertice_start;
	numvert[polygon_num] = 3;
	++polygon_num;
	
	polygon[polygon_num][0] = 3 + vertice_start;
	polygon[polygon_num][1] = 4 + vertice_start;
	polygon[polygon_num][2] = 5 + vertice_start;
	numvert[polygon_num] = 3;
	++polygon_num;
	
	polygon[polygon_num][0] = 0 + vertice_start;
	polygon[polygon_num][1] = 1 + vertice_start;
	polygon[polygon_num][2] = 4 + vertice_start;
	polygon[polygon_num][3] = 3 + vertice_start;
	numvert[polygon_num] = 4;
	++polygon_num;
	
	polygon[polygon_num][0] = 1 + vertice_start;
	polygon[polygon_num][1] = 2 + vertice_start;
	polygon[polygon_num][2] = 5 + vertice_start;
	polygon[polygon_num][3] = 4 + vertice_start;
	numvert[polygon_num] = 4;
	++polygon_num;
	
	polygon[polygon_num][0] = 0 + vertice_start;
	polygon[polygon_num][1] = 2 + vertice_start;
	polygon[polygon_num][2] = 5 + vertice_start;
	polygon[polygon_num][3] = 3 + vertice_start;
	numvert[polygon_num] = 4;
	++polygon_num;
	}
	
	if(zs>=ze)
	{
	polygon[polygon_num][0] = 0 + vertice_start;
	polygon[polygon_num][1] = 1 + vertice_start;
	polygon[polygon_num][2] = 2 + vertice_start;
	numvert[polygon_num] = 3;
	++polygon_num;
	
	polygon[polygon_num][0] = 3 + vertice_start;
	polygon[polygon_num][1] = 4 + vertice_start;
	polygon[polygon_num][2] = 5 + vertice_start;
	numvert[polygon_num] = 3;
	++polygon_num;
	
	polygon[polygon_num][0] = 0 + vertice_start;
	polygon[polygon_num][1] = 1 + vertice_start;
	polygon[polygon_num][2] = 4 + vertice_start;
	polygon[polygon_num][3] = 3 + vertice_start;
	numvert[polygon_num] = 4;
	++polygon_num;
	
	polygon[polygon_num][0] = 0 + vertice_start;
	polygon[polygon_num][1] = 2 + vertice_start;
	polygon[polygon_num][2] = 5 + vertice_start;
	polygon[polygon_num][3] = 3 + vertice_start;
	numvert[polygon_num] = 4;
	++polygon_num;
	
	polygon[polygon_num][0] = 1 + vertice_start;
	polygon[polygon_num][1] = 2 + vertice_start;
	polygon[polygon_num][2] = 5 + vertice_start;
	polygon[polygon_num][3] = 4 + vertice_start;
	numvert[polygon_num] = 4;
	++polygon_num;
	}  
}

void print_porous::wedge_y(lexer *p, fdm *a, ghostcell *pgc,int rank)
{	
	double xs,ys,zs,xe,ye,ze;
	int vertice_start=vertice_num;
	

	xs = p->B282_xs[rank];
    xe = p->B282_xe[rank];

    ys = p->B282_ys[rank];
    ye = p->B282_ye[rank];

    zs = p->B282_zs[rank];
    ze = p->B282_ze[rank];  

	if(zs<ze)
	{
	// Vertices	
		vertice[vertice_num][0] = xs;
		vertice[vertice_num][1] = ys;
		vertice[vertice_num][2] = zs;
		++vertice_num;
		
		vertice[vertice_num][0] = xs;
		vertice[vertice_num][1] = ye;
		vertice[vertice_num][2] = zs;
		++vertice_num;
		
		vertice[vertice_num][0] = xs;
		vertice[vertice_num][1] = ye;
		vertice[vertice_num][2] = ze;
		++vertice_num;
		
		vertice[vertice_num][0] = xe;
		vertice[vertice_num][1] = ys;
		vertice[vertice_num][2] = zs;
		++vertice_num;
		
		vertice[vertice_num][0] = xe;
		vertice[vertice_num][1] = ye;
		vertice[vertice_num][2] = zs;
		++vertice_num;
		
		vertice[vertice_num][0] = xe;
		vertice[vertice_num][1] = ye;
		vertice[vertice_num][2] = ze;
		++vertice_num;
		

	}
	
	if(zs>=ze)
	{
	// Vertices	
		vertice[vertice_num][0] = xs;
		vertice[vertice_num][1] = ys;
		vertice[vertice_num][2] = ze;
		++vertice_num;
		
		vertice[vertice_num][0] = xe;
		vertice[vertice_num][1] = ys;
		vertice[vertice_num][2] = ze;
		++vertice_num;
		
		vertice[vertice_num][0] = xs;
		vertice[vertice_num][1] = ys;
		vertice[vertice_num][2] = zs;
		++vertice_num;
		
		vertice[vertice_num][0] = xs;
		vertice[vertice_num][1] = ye;
		vertice[vertice_num][2] = ze;
		++vertice_num;
		
		vertice[vertice_num][0] = xe;
		vertice[vertice_num][1] = ye;
		vertice[vertice_num][2] = ze;
		++vertice_num;
		
		vertice[vertice_num][0] = xs;
		vertice[vertice_num][1] = ye;
		vertice[vertice_num][2] = zs;
		++vertice_num;
	}
	
		
// Polygon
	if(zs<ze)
	{
	polygon[polygon_num][0] = 0 + vertice_start;
	polygon[polygon_num][1] = 1 + vertice_start;
	polygon[polygon_num][2] = 2 + vertice_start;
	numvert[polygon_num] = 3;
	++polygon_num;
	
	polygon[polygon_num][0] = 3 + vertice_start;
	polygon[polygon_num][1] = 4 + vertice_start;
	polygon[polygon_num][2] = 5 + vertice_start;
	numvert[polygon_num] = 3;
	++polygon_num;
	
	polygon[polygon_num][0] = 0 + vertice_start;
	polygon[polygon_num][1] = 1 + vertice_start;
	polygon[polygon_num][2] = 4 + vertice_start;
	polygon[polygon_num][3] = 3 + vertice_start;
	numvert[polygon_num] = 4;
	++polygon_num;
	
	polygon[polygon_num][0] = 1 + vertice_start;
	polygon[polygon_num][1] = 2 + vertice_start;
	polygon[polygon_num][2] = 5 + vertice_start;
	polygon[polygon_num][3] = 4 + vertice_start;
	numvert[polygon_num] = 4;
	++polygon_num;
	
	polygon[polygon_num][0] = 0 + vertice_start;
	polygon[polygon_num][1] = 2 + vertice_start;
	polygon[polygon_num][2] = 5 + vertice_start;
	polygon[polygon_num][3] = 3 + vertice_start;
	numvert[polygon_num] = 4;
	++polygon_num;
	}
	
	if(zs>=ze)
	{
	polygon[polygon_num][0] = 0 + vertice_start;
	polygon[polygon_num][1] = 1 + vertice_start;
	polygon[polygon_num][2] = 2 + vertice_start;
	numvert[polygon_num] = 3;
	++polygon_num;
	
	polygon[polygon_num][0] = 3 + vertice_start;
	polygon[polygon_num][1] = 4 + vertice_start;
	polygon[polygon_num][2] = 5 + vertice_start;
	numvert[polygon_num] = 3;
	++polygon_num;
	
	polygon[polygon_num][0] = 0 + vertice_start;
	polygon[polygon_num][1] = 1 + vertice_start;
	polygon[polygon_num][2] = 4 + vertice_start;
	polygon[polygon_num][3] = 3 + vertice_start;
	numvert[polygon_num] = 4;
	++polygon_num;
	
	polygon[polygon_num][0] = 0 + vertice_start;
	polygon[polygon_num][1] = 2 + vertice_start;
	polygon[polygon_num][2] = 5 + vertice_start;
	polygon[polygon_num][3] = 3 + vertice_start;
	numvert[polygon_num] = 4;
	++polygon_num;
	
	polygon[polygon_num][0] = 1 + vertice_start;
	polygon[polygon_num][1] = 2 + vertice_start;
	polygon[polygon_num][2] = 5 + vertice_start;
	polygon[polygon_num][3] = 4 + vertice_start;
	numvert[polygon_num] = 4;
	++polygon_num;
	}  
}


void print_porous::plate_x(lexer *p, fdm *a, ghostcell *pgc,int rank)
{	
	double xs,ys,zs,xe,ye,ze,d;
	int vertice_start=vertice_num;
	

	xs = p->B291_xs[rank];
    xe = p->B291_xe[rank];

    ys = p->B291_ys[rank];
    ye = p->B291_ye[rank];

    zs = p->B291_zs[rank];
    ze = p->B291_ze[rank];

    d = p->B291_d[rank];  


	// Vertices	
		vertice[vertice_num][0] = xs;
		vertice[vertice_num][1] = ys;
		vertice[vertice_num][2] = zs;
		++vertice_num;
		
		vertice[vertice_num][0] = xe;
		vertice[vertice_num][1] = ys;
		vertice[vertice_num][2] = ze;
		++vertice_num;
		
		vertice[vertice_num][0] = xe;
		vertice[vertice_num][1] = ye;
		vertice[vertice_num][2] = ze;
		++vertice_num;
		
		vertice[vertice_num][0] = xs;
		vertice[vertice_num][1] = ye;
		vertice[vertice_num][2] = zs;
		++vertice_num;
		
		vertice[vertice_num][0] = xs;
		vertice[vertice_num][1] = ys;
		vertice[vertice_num][2] = zs+d;
		++vertice_num;
		
		vertice[vertice_num][0] = xe;
		vertice[vertice_num][1] = ys;
		vertice[vertice_num][2] = ze+d;
		++vertice_num;
		
		vertice[vertice_num][0] = xe;
		vertice[vertice_num][1] = ye;
		vertice[vertice_num][2] = ze+d;
		++vertice_num;
		
		vertice[vertice_num][0] = xs;
		vertice[vertice_num][1] = ye;
		vertice[vertice_num][2] = zs+d;
		++vertice_num;
	
		
// Polygon
    // bottom
	polygon[polygon_num][0] = 0 + vertice_start;
	polygon[polygon_num][1] = 1 + vertice_start;
	polygon[polygon_num][2] = 2 + vertice_start;
	polygon[polygon_num][3] = 3 + vertice_start;
	numvert[polygon_num] = 4;
	++polygon_num;
	
    // top
	polygon[polygon_num][0] = 4 + vertice_start;
	polygon[polygon_num][1] = 5 + vertice_start;
	polygon[polygon_num][2] = 6 + vertice_start;
	polygon[polygon_num][3] = 7 + vertice_start;
	numvert[polygon_num] = 4;
	++polygon_num;
	
    // front
	polygon[polygon_num][0] = 0 + vertice_start;
	polygon[polygon_num][1] = 4 + vertice_start;
	polygon[polygon_num][2] = 7 + vertice_start;
	polygon[polygon_num][3] = 3 + vertice_start;
	numvert[polygon_num] = 4;
	++polygon_num;
    
    //back
    polygon[polygon_num][0] = 1 + vertice_start;
	polygon[polygon_num][1] = 2 + vertice_start;
	polygon[polygon_num][2] = 6 + vertice_start;
	polygon[polygon_num][3] = 5 + vertice_start;
	numvert[polygon_num] = 4;
	++polygon_num;

    // left
    polygon[polygon_num][0] = 2 + vertice_start;
	polygon[polygon_num][1] = 3 + vertice_start;
	polygon[polygon_num][2] = 7 + vertice_start;
	polygon[polygon_num][3] = 6 + vertice_start;
	numvert[polygon_num] = 4;
	++polygon_num;
    
    // right
    polygon[polygon_num][0] = 0 + vertice_start;
	polygon[polygon_num][1] = 1 + vertice_start;
	polygon[polygon_num][2] = 5 + vertice_start;
	polygon[polygon_num][3] = 4 + vertice_start;
	numvert[polygon_num] = 4;
	++polygon_num;
    
}


void print_porous::box_veg(lexer *p, fdm *a, ghostcell *pgc,int rank)
{
	double xs,ys,zs,xe,ye,ze;
	int vertice_start=vertice_num;
	

	xs = p->B310_xs[rank];
    xe = p->B310_xe[rank];

    ys = p->B310_ys[rank];
    ye = p->B310_ye[rank];

    zs = p->B310_zs[rank];
    ze = p->B310_ze[rank];    
	

	// Vertices	
	vertice[vertice_num][0] = xs;
	vertice[vertice_num][1] = ys;
	vertice[vertice_num][2] = zs;
	++vertice_num;
	
	vertice[vertice_num][0] = xe;
	vertice[vertice_num][1] = ys;
	vertice[vertice_num][2] = zs;
	++vertice_num;
	
	vertice[vertice_num][0] = xe;
	vertice[vertice_num][1] = ye;
	vertice[vertice_num][2] = zs;
	++vertice_num;
	
	vertice[vertice_num][0] = xs;
	vertice[vertice_num][1] = ye;
	vertice[vertice_num][2] = zs;
	++vertice_num;
	
	vertice[vertice_num][0] = xs;
	vertice[vertice_num][1] = ys;
	vertice[vertice_num][2] = ze;
	++vertice_num;
	
	vertice[vertice_num][0] = xe;
	vertice[vertice_num][1] = ys;
	vertice[vertice_num][2] = ze;
	++vertice_num;
	
	vertice[vertice_num][0] = xe;
	vertice[vertice_num][1] = ye;
	vertice[vertice_num][2] = ze;
	++vertice_num;
	
	vertice[vertice_num][0] = xs;
	vertice[vertice_num][1] = ye;
	vertice[vertice_num][2] = ze;
	++vertice_num;
	
	
		
// Polygon

	polygon[polygon_num][0] = 0 + vertice_start;
	polygon[polygon_num][1] = 1 + vertice_start;
	polygon[polygon_num][2] = 2 + vertice_start;
	polygon[polygon_num][3] = 3 + vertice_start;
	numvert[polygon_num] = 4;
	++polygon_num;
	
	polygon[polygon_num][0] = 0 + vertice_start;
	polygon[polygon_num][1] = 1 + vertice_start;
	polygon[polygon_num][2] = 5 + vertice_start;
	polygon[polygon_num][3] = 4 + vertice_start;
	numvert[polygon_num] = 4;
	++polygon_num;
	
	polygon[polygon_num][0] = 1 + vertice_start;
	polygon[polygon_num][1] = 2 + vertice_start;
	polygon[polygon_num][2] = 6 + vertice_start;
	polygon[polygon_num][3] = 5 + vertice_start;
	numvert[polygon_num] = 4;
	++polygon_num;
	
	polygon[polygon_num][0] = 2 + vertice_start;
	polygon[polygon_num][1] = 3 + vertice_start;
	polygon[polygon_num][2] = 7 + vertice_start;
	polygon[polygon_num][3] = 6 + vertice_start;
	numvert[polygon_num] = 4;
	++polygon_num;
	
	polygon[polygon_num][0] = 3 + vertice_start;
	polygon[polygon_num][1] = 0 + vertice_start;
	polygon[polygon_num][2] = 4 + vertice_start;
	polygon[polygon_num][3] = 7 + vertice_start;
	numvert[polygon_num] = 4;
	++polygon_num;
	
	polygon[polygon_num][0] = 4 + vertice_start;
	polygon[polygon_num][1] = 5 + vertice_start;
	polygon[polygon_num][2] = 6 + vertice_start;
	polygon[polygon_num][3] = 7 + vertice_start;
	numvert[polygon_num] = 4;
	++polygon_num;
}

void print_porous::wedge_x_veg(lexer *p, fdm *a, ghostcell *pgc,int rank)
{	
	double xs,ys,zs,xe,ye,ze;
	int vertice_start=vertice_num;
	

	xs = p->B321_xs[rank];
    xe = p->B321_xe[rank];

    ys = p->B321_ys[rank];
    ye = p->B321_ye[rank];

    zs = p->B321_zs[rank];
    ze = p->B321_ze[rank];  

	if(zs<ze)
	{
	// Vertices	
		vertice[vertice_num][0] = xs;
		vertice[vertice_num][1] = ys;
		vertice[vertice_num][2] = zs;
		++vertice_num;
		
		vertice[vertice_num][0] = xe;
		vertice[vertice_num][1] = ys;
		vertice[vertice_num][2] = zs;
		++vertice_num;
		
		vertice[vertice_num][0] = xe;
		vertice[vertice_num][1] = ys;
		vertice[vertice_num][2] = ze;
		++vertice_num;
		
		vertice[vertice_num][0] = xs;
		vertice[vertice_num][1] = ye;
		vertice[vertice_num][2] = zs;
		++vertice_num;
		
		vertice[vertice_num][0] = xe;
		vertice[vertice_num][1] = ye;
		vertice[vertice_num][2] = zs;
		++vertice_num;
		
		vertice[vertice_num][0] = xe;
		vertice[vertice_num][1] = ye;
		vertice[vertice_num][2] = ze;
		++vertice_num;
		

	}
	
	if(zs>=ze)
	{
	// Vertices	
		vertice[vertice_num][0] = xs;
		vertice[vertice_num][1] = ys;
		vertice[vertice_num][2] = ze;
		++vertice_num;
		
		vertice[vertice_num][0] = xe;
		vertice[vertice_num][1] = ys;
		vertice[vertice_num][2] = ze;
		++vertice_num;
		
		vertice[vertice_num][0] = xs;
		vertice[vertice_num][1] = ys;
		vertice[vertice_num][2] = zs;
		++vertice_num;
		
		vertice[vertice_num][0] = xs;
		vertice[vertice_num][1] = ye;
		vertice[vertice_num][2] = ze;
		++vertice_num;
		
		vertice[vertice_num][0] = xe;
		vertice[vertice_num][1] = ye;
		vertice[vertice_num][2] = ze;
		++vertice_num;
		
		vertice[vertice_num][0] = xs;
		vertice[vertice_num][1] = ye;
		vertice[vertice_num][2] = zs;
		++vertice_num;
	}
	
		
// Polygon
	if(zs<ze)
	{
	polygon[polygon_num][0] = 0 + vertice_start;
	polygon[polygon_num][1] = 1 + vertice_start;
	polygon[polygon_num][2] = 2 + vertice_start;
	numvert[polygon_num] = 3;
	++polygon_num;
	
	polygon[polygon_num][0] = 3 + vertice_start;
	polygon[polygon_num][1] = 4 + vertice_start;
	polygon[polygon_num][2] = 5 + vertice_start;
	numvert[polygon_num] = 3;
	++polygon_num;
	
	polygon[polygon_num][0] = 0 + vertice_start;
	polygon[polygon_num][1] = 1 + vertice_start;
	polygon[polygon_num][2] = 4 + vertice_start;
	polygon[polygon_num][3] = 3 + vertice_start;
	numvert[polygon_num] = 4;
	++polygon_num;
	
	polygon[polygon_num][0] = 1 + vertice_start;
	polygon[polygon_num][1] = 2 + vertice_start;
	polygon[polygon_num][2] = 5 + vertice_start;
	polygon[polygon_num][3] = 4 + vertice_start;
	numvert[polygon_num] = 4;
	++polygon_num;
	
	polygon[polygon_num][0] = 0 + vertice_start;
	polygon[polygon_num][1] = 2 + vertice_start;
	polygon[polygon_num][2] = 5 + vertice_start;
	polygon[polygon_num][3] = 3 + vertice_start;
	numvert[polygon_num] = 4;
	++polygon_num;
	}
	
	if(zs>=ze)
	{
	polygon[polygon_num][0] = 0 + vertice_start;
	polygon[polygon_num][1] = 1 + vertice_start;
	polygon[polygon_num][2] = 2 + vertice_start;
	numvert[polygon_num] = 3;
	++polygon_num;
	
	polygon[polygon_num][0] = 3 + vertice_start;
	polygon[polygon_num][1] = 4 + vertice_start;
	polygon[polygon_num][2] = 5 + vertice_start;
	numvert[polygon_num] = 3;
	++polygon_num;
	
	polygon[polygon_num][0] = 0 + vertice_start;
	polygon[polygon_num][1] = 1 + vertice_start;
	polygon[polygon_num][2] = 4 + vertice_start;
	polygon[polygon_num][3] = 3 + vertice_start;
	numvert[polygon_num] = 4;
	++polygon_num;
	
	polygon[polygon_num][0] = 0 + vertice_start;
	polygon[polygon_num][1] = 2 + vertice_start;
	polygon[polygon_num][2] = 5 + vertice_start;
	polygon[polygon_num][3] = 3 + vertice_start;
	numvert[polygon_num] = 4;
	++polygon_num;
	
	polygon[polygon_num][0] = 1 + vertice_start;
	polygon[polygon_num][1] = 2 + vertice_start;
	polygon[polygon_num][2] = 5 + vertice_start;
	polygon[polygon_num][3] = 4 + vertice_start;
	numvert[polygon_num] = 4;
	++polygon_num;
	}  
}

void print_porous::wedge_y_veg(lexer *p, fdm *a, ghostcell *pgc,int rank)
{	
	double xs,ys,zs,xe,ye,ze;
	int vertice_start=vertice_num;
	

	xs = p->B322_xs[rank];
    xe = p->B322_xe[rank];

    ys = p->B322_ys[rank];
    ye = p->B322_ye[rank];

    zs = p->B322_zs[rank];
    ze = p->B322_ze[rank];  

	if(zs<ze)
	{
	// Vertices	
		vertice[vertice_num][0] = xs;
		vertice[vertice_num][1] = ys;
		vertice[vertice_num][2] = zs;
		++vertice_num;
		
		vertice[vertice_num][0] = xs;
		vertice[vertice_num][1] = ye;
		vertice[vertice_num][2] = zs;
		++vertice_num;
		
		vertice[vertice_num][0] = xs;
		vertice[vertice_num][1] = ye;
		vertice[vertice_num][2] = ze;
		++vertice_num;
		
		vertice[vertice_num][0] = xe;
		vertice[vertice_num][1] = ys;
		vertice[vertice_num][2] = zs;
		++vertice_num;
		
		vertice[vertice_num][0] = xe;
		vertice[vertice_num][1] = ye;
		vertice[vertice_num][2] = zs;
		++vertice_num;
		
		vertice[vertice_num][0] = xe;
		vertice[vertice_num][1] = ye;
		vertice[vertice_num][2] = ze;
		++vertice_num;
		

	}
	
	if(zs>=ze)
	{
	// Vertices	
		vertice[vertice_num][0] = xs;
		vertice[vertice_num][1] = ys;
		vertice[vertice_num][2] = ze;
		++vertice_num;
		
		vertice[vertice_num][0] = xe;
		vertice[vertice_num][1] = ys;
		vertice[vertice_num][2] = ze;
		++vertice_num;
		
		vertice[vertice_num][0] = xs;
		vertice[vertice_num][1] = ys;
		vertice[vertice_num][2] = zs;
		++vertice_num;
		
		vertice[vertice_num][0] = xs;
		vertice[vertice_num][1] = ye;
		vertice[vertice_num][2] = ze;
		++vertice_num;
		
		vertice[vertice_num][0] = xe;
		vertice[vertice_num][1] = ye;
		vertice[vertice_num][2] = ze;
		++vertice_num;
		
		vertice[vertice_num][0] = xs;
		vertice[vertice_num][1] = ye;
		vertice[vertice_num][2] = zs;
		++vertice_num;
	}
	
		
// Polygon
	if(zs<ze)
	{
	polygon[polygon_num][0] = 0 + vertice_start;
	polygon[polygon_num][1] = 1 + vertice_start;
	polygon[polygon_num][2] = 2 + vertice_start;
	numvert[polygon_num] = 3;
	++polygon_num;
	
	polygon[polygon_num][0] = 3 + vertice_start;
	polygon[polygon_num][1] = 4 + vertice_start;
	polygon[polygon_num][2] = 5 + vertice_start;
	numvert[polygon_num] = 3;
	++polygon_num;
	
	polygon[polygon_num][0] = 0 + vertice_start;
	polygon[polygon_num][1] = 1 + vertice_start;
	polygon[polygon_num][2] = 4 + vertice_start;
	polygon[polygon_num][3] = 3 + vertice_start;
	numvert[polygon_num] = 4;
	++polygon_num;
	
	polygon[polygon_num][0] = 1 + vertice_start;
	polygon[polygon_num][1] = 2 + vertice_start;
	polygon[polygon_num][2] = 5 + vertice_start;
	polygon[polygon_num][3] = 4 + vertice_start;
	numvert[polygon_num] = 4;
	++polygon_num;
	
	polygon[polygon_num][0] = 0 + vertice_start;
	polygon[polygon_num][1] = 2 + vertice_start;
	polygon[polygon_num][2] = 5 + vertice_start;
	polygon[polygon_num][3] = 3 + vertice_start;
	numvert[polygon_num] = 4;
	++polygon_num;
	}
	
	if(zs>=ze)
	{
	polygon[polygon_num][0] = 0 + vertice_start;
	polygon[polygon_num][1] = 1 + vertice_start;
	polygon[polygon_num][2] = 2 + vertice_start;
	numvert[polygon_num] = 3;
	++polygon_num;
	
	polygon[polygon_num][0] = 3 + vertice_start;
	polygon[polygon_num][1] = 4 + vertice_start;
	polygon[polygon_num][2] = 5 + vertice_start;
	numvert[polygon_num] = 3;
	++polygon_num;
	
	polygon[polygon_num][0] = 0 + vertice_start;
	polygon[polygon_num][1] = 1 + vertice_start;
	polygon[polygon_num][2] = 4 + vertice_start;
	polygon[polygon_num][3] = 3 + vertice_start;
	numvert[polygon_num] = 4;
	++polygon_num;
	
	polygon[polygon_num][0] = 0 + vertice_start;
	polygon[polygon_num][1] = 2 + vertice_start;
	polygon[polygon_num][2] = 5 + vertice_start;
	polygon[polygon_num][3] = 3 + vertice_start;
	numvert[polygon_num] = 4;
	++polygon_num;
	
	polygon[polygon_num][0] = 1 + vertice_start;
	polygon[polygon_num][1] = 2 + vertice_start;
	polygon[polygon_num][2] = 5 + vertice_start;
	polygon[polygon_num][3] = 4 + vertice_start;
	numvert[polygon_num] = 4;
	++polygon_num;
	}  
}
