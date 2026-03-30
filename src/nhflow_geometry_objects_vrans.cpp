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

#include"nhflow_geometry.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

void nhflow_geometry::objects_create_vrans(lexer *p, ghostcell *pgc)
{
    int qn;

    objects_allocate_vrans(p,pgc);
	
    entity_count=0;
	
	for(qn=0;qn<p->B210;++qn)
    {
        box(p,pgc,qn);
        ++entity_count;
    }
    
    for(qn=0;qn<p->B212;++qn)
    {
        cylinder_y(p,pgc,qn);
        ++entity_count;
    }
    
	for(qn=0;qn<p->B213;++qn)
    {
        cylinder_z(p,pgc,qn);
        ++entity_count;
    }
    
    for(qn=0;qn<p->B214;++qn)
    {
        jacketmember(p,pgc,qn);
        ++entity_count;
    }
    
    for(qn=0;qn<p->B215;++qn)
    {
        sphere(p,pgc,qn);
        ++entity_count;
    }
    
    for(qn=0;qn<p->B217;++qn)
    {
        wedge_x(p,pgc,qn);
        ++entity_count;
    }
    
    for(qn=0;qn<p->B218;++qn)
    {
        wedge_y(p,pgc,qn);
        ++entity_count;
    }
    
    for(qn=0;qn<p->B219;++qn)
    {
        wedge_z(p,pgc,qn);
        ++entity_count;
    }
    
    if(p->B230==1)
    {
        read_stl(p,pgc);
		++entity_count;
    }
    
    if(p->mpirank==0)
    print_vtp(p,1);

    if(p->mpirank==0)
	cout<<"Surface triangles: "<<tricount<<endl;
}

void nhflow_geometry::objects_allocate_vrans(lexer *p, ghostcell *pgc)
{
    double U,ds,phi,r,snum,trisum;
    
    entity_sum = p->B210 + p->B212 + p->B213 + p->B214 + p->B215 + p->B217 + p->B218 + p->B219;
	tricount=0;
    trisum=0;
    
    // box
    trisum+=12*p->B210;
    
    // cylinder_y
    for(n=0; n<p->B212;++n)
	{
	r = p->B212_r[n];
	U = 2.0*PI*r;
	ds = 0.5*(DSM);
    snum = MAX(int(U/ds), 40);
	trisum+=6*snum;
	}
    
    // cylinder_z
    for(n=0; n<p->B213;++n)
	{
	r = p->B213_r[n];
	U = 2.0*PI*r;
	ds = 0.5*(DSM);
    snum = MAX(int(U/ds), 40);
	trisum+=6*snum;
	}
    
    // cylinder_member
    for(n=0; n<p->B214;++n)
	{
	r = MAX(p->B214_r1[n],p->B214_r2[n]);
	U = 2.0*PI*r;
	ds = 0.5*(DSM);
    snum = MAX(int(U/ds), 40);
	trisum+=6*snum;
	}
    
    // sphere
    for(n=0; n<p->B215;++n)
    {
	r = p->B215_r[n];
	U = 2.0*PI*r;
	ds = 0.85*(DSM);
    snum = MAX(int(U/ds), 40);
    trisum+=snum*snum*2;
    }
    
    // wedge
    trisum+=8*p->B217;
    
    // wedge
    trisum+=8*p->B218;
    
    // wedge
    trisum+=8*p->B219;

    // STL
    if(p->B230==1)
    entity_sum+=1;

    p->Darray(tri_x,trisum,3);
	p->Darray(tri_y,trisum,3);
	p->Darray(tri_z,trisum,3);
    p->Darray(tri_x0,trisum,3);
	p->Darray(tri_y0,trisum,3);
	p->Darray(tri_z0,trisum,3);   
    
	p->Iarray(tstart,entity_sum);
	p->Iarray(tend,entity_sum);
}
