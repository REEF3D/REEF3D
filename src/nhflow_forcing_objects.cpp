/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"nhflow_forcing.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

void nhflow_forcing::objects_create(lexer *p, ghostcell *pgc)
{
    int qn;

    objects_allocate(p,pgc);
	
    entity_count=0;
	
	for(qn=0;qn<p->A581;++qn)
    {
        box(p,pgc,qn);
        ++entity_count;
    }
    
    for(qn=0;qn<p->A583;++qn)
    {
        cylinder_y(p,pgc,qn);
        ++entity_count;
    }
    
	for(qn=0;qn<p->A584;++qn)
    {
        cylinder_z(p,pgc,qn);
        ++entity_count;
    }
    
    for(qn=0;qn<p->A585;++qn)
    {
        jacketmember(p,pgc,qn);
        ++entity_count;
    }
    
    for(qn=0;qn<p->A586;++qn)
    {
        sphere(p,pgc,qn);
        ++entity_count;
    }
    
    for(qn=0;qn<p->A587;++qn)
    {
        wedge_x(p,pgc,qn);
        ++entity_count;
    }
    
    for(qn=0;qn<p->A588;++qn)
    {
        wedge_y(p,pgc,qn);
        ++entity_count;
    }
    
    for(qn=0;qn<p->A589;++qn)
    {
        wedge_z(p,pgc,qn);
        ++entity_count;
    }
    
    
    if(p->A590==1)
    {
        read_stl(p,pgc);
		++entity_count;
    }

    if(p->mpirank==0)
	cout<<"Surface triangles: "<<tricount<<endl;
}

void nhflow_forcing::objects_allocate(lexer *p, ghostcell *pgc)
{
    double U,ds,phi,r,snum,trisum;
    
    entity_sum = p->A581 + p->A583 + p->A584 + p->A585 + p->A586 + p->A587 + p->A588 + p->A589;
	tricount=0;
    trisum=0;
    
    // box
    trisum+=12*p->A581;
    
    // cylinder_y
    for(n=0; n<p->A583;++n)
	{
	r = p->A583_r[n];
	U = 2.0*PI*r;
	ds = 0.75*(U*p->DXM);
	snum = int(U/ds);
	trisum+=6*snum;
	}
    
    // cylinder_z
    for(n=0; n<p->A584;++n)
	{
	r = p->A584_r[n];
	U = 2.0*PI*r;
	ds = 0.75*(U*p->DXM);
	snum = int(U/ds);
	trisum+=6*snum;
	}
    
    // cylinder_member
    for(n=0; n<p->A585;++n)
	{
	r = MAX(p->A585_r1[n],p->A585_r2[n]);
	U = 2.0*PI*r;
	ds = 0.75*(U*p->DXM);
	snum = int(U/ds);
	trisum+=6*snum;
	}
    
    // sphere
    for(n=0; n<p->A586;++n)
    {
	r = p->A586_r[n];
	U = 2.0*PI*r;
	ds = 0.75*(U*p->DXM);
	snum = int(U/ds);
    trisum+=snum*snum*2;
    }
    
    // wedge
    trisum+=8*p->A587;
    
    // wedge
    trisum+=8*p->A588;
    
    // wedge
    trisum+=8*p->A589;

    // STL
    if(p->A590==1)
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
