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

#include"6DOF_gc.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void sixdof_gc::objects(lexer *p, fdm *a, ghostcell *pgc)
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
		cout<<"Multiple floating bodies are not supported yet."<<endl<<endl;
		pgc->final();
		exit(0);
	}

    if(p->mpirank==0)
	cout<<"Surface triangles: "<<tricount<<endl;
	
	// Refine triangles
	//geometry_refinement(p,pgc);	

    if(p->mpirank==0)
	cout<<"Refined surface triangles: "<<tricount<<endl;
}


void sixdof_gc::objects_allocate(lexer *p, fdm *a, ghostcell *pgc)
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
	ds = 0.25*(U*p->DXM);
	snum = int(U/ds);
	trisum+=5*(snum+1)*p->X131;
    
    // cylinder_y
    r=p->X132_rad;
	U = 2.0 * PI * r;
	ds = 0.25*(U*p->DXM);
	snum = int(U/ds);
	trisum+=5*(snum+1)*p->X132;
    
    // cylinder_y
    r=p->X133_rad;
	U = 2.0 * PI * r;
	ds = 0.25*(U*p->DXM);
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
    p->Darray(tri_xn,trisum,3);
	p->Darray(tri_yn,trisum,3);
	p->Darray(tri_zn,trisum,3);    	
    
	p->Iarray(tstart,entity_sum);
	p->Iarray(tend,entity_sum);
}


void sixdof_gc::geometry_refinement(lexer *p, ghostcell *pgc)
{
	double x0,x1,x2,y0,y1,y2,z0,z1,z2;
	double x01,x02,x12,y01,y02,y12,z01,z02,z12;
	double at,bt,ct,st;
	double nx_old,ny_old,nz_old;	
	
	tri_x_r.reserve(4*tricount);
	tri_y_r.reserve(4*tricount);
	tri_z_r.reserve(4*tricount);	
	
	tri_x_r.resize(tricount,vector<double>(3,0.0));
	tri_y_r.resize(tricount,vector<double>(3,0.0));
	tri_z_r.resize(tricount,vector<double>(3,0.0));
	
	
	for (int i = 0; i < tricount; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			tri_x_r[i][j] = tri_x[i][j];
			tri_y_r[i][j] = tri_y[i][j];
			tri_z_r[i][j] = tri_z[i][j];
		}
	}
	
	
	double critL = p->DXM*1.0;
	
	for (int n = 0; n < tri_x_r.size(); n++)
	{
		x0 = tri_x_r[n][0];
		x1 = tri_x_r[n][1];
		x2 = tri_x_r[n][2];
		
		y0 = tri_y_r[n][0];
		y1 = tri_y_r[n][1];
		y2 = tri_y_r[n][2];
		
		z0 = tri_z_r[n][0];
		z1 = tri_z_r[n][1];
		z2 = tri_z_r[n][2];  
           
		at = sqrt(pow(x1-x0,2.0) + pow(y1-y0,2.0) + pow(z1-z0,2.0));
		bt = sqrt(pow(x1-x2,2.0) + pow(y1-y2,2.0) + pow(z1-z2,2.0));
		ct = sqrt(pow(x2-x0,2.0) + pow(y2-y0,2.0) + pow(z2-z0,2.0));   
		   
		
		// Check size of triangle and split into 4 triangles if too big
		
		if ((at + bt + ct)/3.0 > critL)
		{
			// Half points
			
			x01 = x0 + (x1 - x0)/2.0;
			y01 = y0 + (y1 - y0)/2.0;
			z01 = z0 + (z1 - z0)/2.0;

			x02 = x0 + (x2 - x0)/2.0;
			y02 = y0 + (y2 - y0)/2.0;
			z02 = z0 + (z2 - z0)/2.0;			

			x12 = x1 + (x2 - x1)/2.0;
			y12 = y1 + (y2 - y1)/2.0;
			z12 = z1 + (z2 - z1)/2.0;
			
			
			// Old normal vector    
				
			nx_old = (y1 - y0) * (z2 - z0) - (y2 - y0) * (z1 - z0);
			ny_old = (x2 - x0) * (z1 - z0) - (x1 - x0) * (z2 - z0); 
			nz_old = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
			
			
			// Delete old and add new triangles
		
			tri_x_r.erase(tri_x_r.begin() + n); 
			tri_y_r.erase(tri_y_r.begin() + n); 
			tri_z_r.erase(tri_z_r.begin() + n); 
			n--;

			create_triangle(x0,y0,z0,x01,y01,z01,x02,y02,z02,nx_old,ny_old,nz_old);
			create_triangle(x01,y01,z01,x12,y12,z12,x02,y02,z02,nx_old,ny_old,nz_old);
			create_triangle(x01,y01,z01,x1,y1,z1,x12,y12,z12,nx_old,ny_old,nz_old);
			create_triangle(x02,y02,z02,x12,y12,z12,x2,y2,z2,nx_old,ny_old,nz_old);
		}
		
		if (tri_x_r.size() > 20000) break;
	}
	
	
	p->Dresize(tri_x,tricount,tri_x_r.size(),3,3);
	p->Dresize(tri_y,tricount,tri_y_r.size(),3,3);
	p->Dresize(tri_z,tricount,tri_z_r.size(),3,3);
	p->Dresize(tri_xn,tricount,tri_x_r.size(),3,3);
	p->Dresize(tri_yn,tricount,tri_y_r.size(),3,3);
	p->Dresize(tri_zn,tricount,tri_z_r.size(),3,3);
	
	tricount = tri_x_r.size();
	tend[0] = tricount;
	
	for (int i = 0; i < tricount; i++)
	{
		for (int j = 0; j < 3; j++)
		{	
			tri_x[i][j] = tri_x_r[i][j];
			tri_y[i][j] = tri_y_r[i][j];
			tri_z[i][j] = tri_z_r[i][j];
		}
	}
}


void sixdof_gc::create_triangle
(
	double& x0, double& y0, double& z0,
	double& x1, double& y1, double& z1,
	double& x2, double& y2, double& z2,
	const double& nx_old, const double& ny_old, const double& nz_old
)
{
	double nx,ny,nz,temp;
	
	vector<double> tri_x_new(3,0.0);
	vector<double> tri_y_new(3,0.0);
	vector<double> tri_z_new(3,0.0);

	
	// Calculate new normal vector
	
	nx = (y1 - y0) * (z2 - z0) - (y2 - y0) * (z1 - z0);
	ny = (x2 - x0) * (z1 - z0) - (x1 - x0) * (z2 - z0); 
	nz = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);		

	nx = nx > 1.0e-5 ? nx : nx_old;
	ny = ny > 1.0e-5 ? ny : ny_old;
	nz = nz > 1.0e-5 ? nz : nz_old;	
	
	
	// Arrange triangle such that normal vector points outward
	
	if 
	(
		   SIGN(nx) != SIGN(nx_old) 
		|| SIGN(ny) != SIGN(ny_old) 
		|| SIGN(nz) != SIGN(nz_old)
	)
	{
		tri_x_new[0] = x2;
		tri_x_new[1] = x1;
		tri_x_new[2] = x0;

		tri_y_new[0] = y2;
		tri_y_new[1] = y1;
		tri_y_new[2] = y0;

		tri_z_new[0] = z2;
		tri_z_new[1] = z1;
		tri_z_new[2] = z0;				
	}
	else
	{	
		tri_x_new[0] = x0;
		tri_x_new[1] = x1;
		tri_x_new[2] = x2;

		tri_y_new[0] = y0;
		tri_y_new[1] = y1;
		tri_y_new[2] = y2;

		tri_z_new[0] = z0;
		tri_z_new[1] = z1;
		tri_z_new[2] = z2;	
	}
	
	
	// Add triangle to list
	
	tri_x_r.push_back(tri_x_new);
	tri_y_r.push_back(tri_y_new);
	tri_z_r.push_back(tri_z_new);
}
