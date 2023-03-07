/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"6DOF_df_object.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void sixdof_df_object::objects_create(lexer *p, fdm *a, ghostcell *pgc)
{
    int qn;

    objects_allocate(p,a,pgc);
	
    entity_count=0;
	
	for(qn=0;qn<p->X110;++qn)
    {
        box(p,a,pgc,qn);
        ++entity_count;
    }
    
    for(qn=0;qn<p->X131;++qn)
    {
        cylinder_x(p,a,pgc,qn);
        ++entity_count;
    }
	
	for(qn=0;qn<p->X132;++qn)
    {
        cylinder_y(p,a,pgc,qn);
        ++entity_count;
    }
	
	for(qn=0;qn<p->X133;++qn)
    {
        cylinder_z(p,a,pgc,qn);
        ++entity_count;
    }
	
	for(qn=0;qn<p->X153;++qn)
    {
        wedge_sym(p,a,pgc,qn);
        ++entity_count;
    }
    
    for(qn=0;qn<p->X163;++qn)
    {
        wedge(p,a,pgc,qn);
        ++entity_count;
    }
    
    for(qn=0;qn<p->X164;++qn)
    {
        hexahedron(p,a,pgc,qn);
        ++entity_count;
    }
    
    if(p->X180==1)
    {
        read_stl(p,a,pgc);
		++entity_count;
    }


	if (entity_count > 1)
	{
		cout<<"Multiple floating bodies are not fully supported yet."<<endl<<endl;
		//pgc->final();
		//exit(0);
	}

    if(p->mpirank==0)
	cout<<"Surface triangles: "<<tricount<<endl;
    
    // Initialise STL geometric parameters
	geometry_stl(p,a,pgc);
    
    // Order Triangles for correct inside/outside orientation
    triangle_switch_ray(p,a,pgc);
	
	// Refine triangles
	geometry_refinement(p,pgc);	

    if(p->mpirank==0)
	cout<<"Refined surface triangles: "<<tricount<<endl;
}

void sixdof_df_object::objects_allocate(lexer *p, fdm *a, ghostcell *pgc)
{
    double U,ds,phi,r,snum,trisum;
    
    entity_sum = p->X110 + p->X131 + p->X132 + p->X133 + p->X153 + p->X163 + p->X164;
	tricount=0;
    trisum=0;
    
    // box
    trisum+=12*p->X110;
    
    // cylinder_x   
    r=p->X131_rad;
	U = 2.0 * PI * r;
	ds = 0.75*(U*p->dx);
	snum = int(U/ds);
	trisum+=5*(snum+1)*p->X131;
    
    // cylinder_y
    r=p->X132_rad;
	U = 2.0 * PI * r;
	ds = 0.75*(U*p->dx);
	snum = int(U/ds);
	trisum+=5*(snum+1)*p->X132;
    
    // cylinder_y
    r=p->X133_rad;
	U = 2.0 * PI * r;
	ds = 0.75*(U*p->dx);
	snum = int(U/ds);
    trisum+=5*(snum+1)*p->X133;
    
    // wedge sym
    trisum+=12*p->X153;
    
    // wedge
    trisum+=8*p->X163;
    
    // hexahedron
    trisum+=12*p->X164;
    
    // STL
    if(p->X180==1)
    entity_sum=1;

    
    p->Darray(tri_x,trisum,3);
	p->Darray(tri_y,trisum,3);
	p->Darray(tri_z,trisum,3);
    p->Darray(tri_x0,trisum,3);
	p->Darray(tri_y0,trisum,3);
	p->Darray(tri_z0,trisum,3);    	
    
	p->Iarray(tstart,entity_sum);
	p->Iarray(tend,entity_sum);
}
