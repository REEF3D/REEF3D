/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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
--------------------------------------------------------------------*/

#include"print_porous.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"


void print_porous::box(lexer *p, fdm *a, ghostcell *pgc,int rank)
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