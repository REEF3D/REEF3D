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
	
	for(qn=0;qn<p->A561;++qn)
    {
        box(p,pgc,qn);
        ++entity_count;
    }
	
	for(qn=0;qn<p->A564;++qn)
    {
        cylinder_z(p,pgc,qn);
        ++entity_count;
    }
	/*
    if(p->X180==1)
    {
        read_stl(p,pgc);
		++entity_count;
    }*/

    if(p->mpirank==0)
	cout<<"Surface triangles: "<<tricount<<endl;
    
    // Initialise STL geometric parameters
	//geometry_stl(p,pgc);
    
    // Order Triangles for correct inside/outside orientation
    /*if(p->A10==6)
    triangle_switch_ray(p,pgc);
	
	// Refine triangles
    if(p->X185>0 && p->X60==1)
	geometry_refinement(p,pgc);	

    if(p->mpirank==0)
	cout<<"Refined surface triangles: "<<tricount<<endl;*/
}

void nhflow_forcing::objects_allocate(lexer *p, ghostcell *pgc)
{
    double U,ds,phi,r,snum,trisum;
    
    entity_sum = p->A561 + p->A564;
	tricount=0;
    trisum=0;
    
    // box
    trisum+=12*p->A561;
    
    // cylinder_z
    r=p->X133_rad;
	U = 2.0 * PI * r;
	ds = 0.75*(U*p->dx);
	snum = int(U/ds);
    trisum+=5*(snum+1)*p->A564;

    // STL
    //if(p->X180==1)
    //entity_sum=1;

    
    p->Darray(tri_x,trisum,3);
	p->Darray(tri_y,trisum,3);
	p->Darray(tri_z,trisum,3);
   	
    
	p->Iarray(tstart,entity_sum);
	p->Iarray(tend,entity_sum);
}
