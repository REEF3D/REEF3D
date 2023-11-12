/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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
#include<sys/stat.h>
#include<sys/types.h>

print_porous::print_porous(lexer* p, fdm *a, ghostcell *pgc)
{
	vertice_num=0;
	polygon_num=0;
	polygon_sum=0;
	
	vertice_alloc = p->B270*8 + p->B281*6 + p->B282*6 + p->B291*8 + p->B310*8 + p->B321*6 + p->B322*6;
	polygon_alloc = p->B270*6 + p->B281*6 + p->B282*6 + p->B291*6 + p->B310*6 + p->B321*6 + p->B322*6;
    
    for(n=0; n<p->B274;++n)
	{
    
    double U,ds;
    int snum;
	U = 2.0 * PI * p->B274_r[n];
	
	ds = 0.2*(U*p->DXM);
	
	snum = int(U/ds);
	
	vertice_alloc+=12*snum;
    polygon_alloc+=4*snum;
	}
    

	p->Darray(vertice,vertice_alloc,3);
	p->Iarray(polygon,polygon_alloc,4);
	p->Iarray(numvert,polygon_alloc);
}

print_porous::~print_porous()
{
}

void print_porous::start(lexer *p, fdm *a, ghostcell *pgc)
{
	objects(p,a,pgc);
	print_vtp(p,a,pgc);
}

void print_porous::objects(lexer *p, fdm *a, ghostcell *pgc)
{
	int qn;
	
	for(qn=0;qn<p->B270;++qn)
	box(p,a,pgc,qn);
    
    for(qn=0;qn<p->B274;++qn)
	cylinder_z(p,a,pgc,qn);
	
	for(qn=0;qn<p->B281;++qn)
	wedge_x(p,a,pgc,qn);
    
    for(qn=0;qn<p->B282;++qn)
	wedge_y(p,a,pgc,qn);
    
    for(qn=0;qn<p->B291;++qn)
	plate_x(p,a,pgc,qn);
    
    for(qn=0;qn<p->B310;++qn)
	box_veg(p,a,pgc,qn);
    
    for(qn=0;qn<p->B321;++qn)
	wedge_x_veg(p,a,pgc,qn);
    
    for(qn=0;qn<p->B322;++qn)
	wedge_y_veg(p,a,pgc,qn);
	
	for(n=0;n<polygon_num;++n)
	polygon_sum+=numvert[n];
}
