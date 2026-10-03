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

#include"6DOF_geometry.h"
#include"lexer.h"
#include"ghostcell.h"
#include"geo_primitive.h"
#include<math.h>

sixdof_geometry::sixdof_geometry(int number) : id(number), tri_x(nullptr), tri_y(nullptr), tri_z(nullptr),
                                               tri_x0(nullptr), tri_y0(nullptr), tri_z0(nullptr),
                                               tricount(0), entity_count(0), entity_sum(0),
                                               tstart(nullptr), tend(nullptr), holes_closed(false), amr_hfac(1.0),
                                               tri_switch(nullptr), tri_switch_local(nullptr),
                                               tricount_local(0), tricount_local_list(nullptr), tricount_local_displ(nullptr),
                                               tricount_switch_total(0), n(0), q(0)
{
}

sixdof_geometry::~sixdof_geometry()
{
}

// ---------------------------------------------------------------------------------------------
// set-up
// ---------------------------------------------------------------------------------------------

void sixdof_geometry::allocate(lexer *p)
{
    double U,ds,r,snum,trisum;
    
    entity_sum = p->X110 + p->X131 + p->X132 + p->X133 + p->X153 + p->X163 + p->X164 + p->X170 + p->X171 + p->X172;
    
	tricount=0;
    entity_count=0;
    trisum=0;
    
    // box
    trisum+=12*p->X110;
    
    // cylinder_x   
    r=p->X131_rad;
	U = 2.0 * PI * r;
	ds = 0.75*p->dx;      // segment length 0.75 dx (was 0.75*U*dx, i.e. 1/(0.75 dx) segments for any radius)
	snum = MAX(int(U/ds),8);
	trisum+=5*(snum+1)*p->X131;
    
    // cylinder_y
    r=p->X132_rad;
	U = 2.0 * PI * r;
	ds = 0.75*p->dx;      // segment length 0.75 dx (was 0.75*U*dx, i.e. 1/(0.75 dx) segments for any radius)
	snum = MAX(int(U/ds),8);
	trisum+=5*(snum+1)*p->X132;
    
    // cylinder_z
    r=p->X133_rad;
	U = 2.0 * PI * r;
	ds = 0.75*p->dx;      // segment length 0.75 dx (was 0.75*U*dx, i.e. 1/(0.75 dx) segments for any radius)
	snum = MAX(int(U/ds),8);
    trisum+=5*(snum+1)*p->X133;
    
    // wedge sym
    trisum+=12*p->X153;
    
    // wedge
    trisum+=8*p->X163;
    
    // hexahedron
    trisum+=12*p->X164;
    
    // piston
    trisum+=12*p->X170;
    
    // double flap
    trisum+=28*p->X172;
    
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

void sixdof_geometry::primitives(lexer *p, ghostcell *pgc)
{
    int qn;
    
	for(qn=0;qn<p->X110;++qn)
    {
        box(p,qn);
        ++entity_count;
    }
	
    for(qn=0;qn<p->X131;++qn)
    {
        cylinder_x(p,qn);
        ++entity_count;
    }
	
	for(qn=0;qn<p->X132;++qn)
    {
        cylinder_y(p,qn);
        ++entity_count;
    }
	
	for(qn=0;qn<p->X133;++qn)
    {
        cylinder_z(p,qn);
        ++entity_count;
    }
	
	for(qn=0;qn<p->X153;++qn)
    {
        wedge_sym(p,qn);
        ++entity_count;
    }
	
    for(qn=0;qn<p->X163;++qn)
    {
        wedge(p,qn);
        ++entity_count;
    }
    
    for(qn=0;qn<p->X164;++qn)
    {
        hexahedron(p,qn);
        ++entity_count;
    }
}

void sixdof_geometry::refine(lexer *p, ghostcell *pgc)
{
    // X 185 1-3: split triangles larger than the grid
    // X 185 4  : adaptive isotropic re-triangulation matching the local grid spacing
    if(p->X185==4)
    geometry_remesh(p,pgc);
    
    else
    geometry_refinement(p,pgc);
}

// ---------------------------------------------------------------------------------------------
// primitives
// ---------------------------------------------------------------------------------------------

void sixdof_geometry::box(lexer *p, int qn)
{
	tstart[entity_count]=tricount;
    geo_primitive::box(tri_x,tri_y,tri_z,tricount,p->X110_xs[qn],p->X110_xe[qn],p->X110_ys[qn],p->X110_ye[qn],p->X110_zs[qn],p->X110_ze[qn]);
	tend[entity_count]=tricount;
}

void sixdof_geometry::cylinder_x(lexer *p, int qn)
{
    // segment length 0.75 dx (was 0.75*U*dx, i.e. 1/(0.75 dx) segments for any radius)
	const int snum = geo_primitive::segments(p->X131_rad,0.75*p->dx,8);
	tstart[entity_count]=tricount;
	geo_primitive::cylinder_x(tri_x,tri_y,tri_z,tricount,p->X131_yc,p->X131_zc,p->X131_xc-0.5*p->X131_h,p->X131_xc+0.5*p->X131_h,p->X131_rad,snum);
	tend[entity_count]=tricount;
}

void sixdof_geometry::cylinder_y(lexer *p, int qn)
{
	const int snum = geo_primitive::segments(p->X132_rad,0.75*p->dx,8);
	tstart[entity_count]=tricount;
	geo_primitive::cylinder_y(tri_x,tri_y,tri_z,tricount,p->X132_xc,p->X132_zc,p->X132_yc-0.5*p->X132_h,p->X132_yc+0.5*p->X132_h,p->X132_rad,snum);
	tend[entity_count]=tricount;
}

void sixdof_geometry::cylinder_z(lexer *p, int qn)
{
	const int snum = geo_primitive::segments(p->X133_rad,0.75*p->dx,8);
	tstart[entity_count]=tricount;
	geo_primitive::cylinder_z(tri_x,tri_y,tri_z,tricount,p->X133_xc,p->X133_yc,p->X133_zc-0.5*p->X133_h,p->X133_zc+0.5*p->X133_h,p->X133_rad,snum);
	tend[entity_count]=tricount;
}

void sixdof_geometry::wedge_sym(lexer *p, int qn)
{
	tstart[entity_count]=tricount;
    geo_primitive::wedge_sym(tri_x,tri_y,tri_z,tricount,p->X153_xs,p->X153_xe,p->X153_ys,p->X153_ye,p->X153_zs,p->X153_ze);
	tend[entity_count]=tricount;
}

void sixdof_geometry::wedge(lexer *p, int qn)
{
    const double v[18] = {p->X163_x1[qn], p->X163_y1[qn], p->X163_z1[qn], p->X163_x2[qn], p->X163_y2[qn], p->X163_z2[qn], p->X163_x3[qn], p->X163_y3[qn], p->X163_z3[qn], p->X163_x4[qn], p->X163_y4[qn], p->X163_z4[qn], p->X163_x5[qn], p->X163_y5[qn], p->X163_z5[qn], p->X163_x6[qn], p->X163_y6[qn], p->X163_z6[qn]};
	tstart[entity_count]=tricount;
    geo_primitive::wedge(tri_x,tri_y,tri_z,tricount,v);
    tend[entity_count]=tricount;
}

void sixdof_geometry::hexahedron(lexer *p, int qn)
{
    const double v[24] = {p->X164_x1[qn], p->X164_y1[qn], p->X164_z1[qn], p->X164_x2[qn], p->X164_y2[qn], p->X164_z2[qn], p->X164_x3[qn], p->X164_y3[qn], p->X164_z3[qn], p->X164_x4[qn], p->X164_y4[qn], p->X164_z4[qn], p->X164_x5[qn], p->X164_y5[qn], p->X164_z5[qn], p->X164_x6[qn], p->X164_y6[qn], p->X164_z6[qn], p->X164_x7[qn], p->X164_y7[qn], p->X164_z7[qn], p->X164_x8[qn], p->X164_y8[qn], p->X164_z8[qn]};
	tstart[entity_count]=tricount;
    geo_primitive::hexahedron(tri_x,tri_y,tri_z,tricount,v);
    tend[entity_count]=tricount;
}

// ---------------------------------------------------------------------------------------------
// pose
// ---------------------------------------------------------------------------------------------

void sixdof_geometry::store_body_frame(const Eigen::Vector3d &c)
{
	for(n=0; n<tricount; ++n)
	{
        for(int q=0; q<3; q++)
        {        
            tri_x0[n][q] = tri_x[n][q] - c(0);
            tri_y0[n][q] = tri_y[n][q] - c(1);
            tri_z0[n][q] = tri_z[n][q] - c(2);
        }
    }
}

void sixdof_geometry::rotate(double phi, double theta, double psi, const Eigen::Vector3d &c)
{
    for (n=0; n<tricount; ++n)
    for (int q=0; q<3; ++q)
    rotation_tri(phi,theta,psi,tri_x[n][q],tri_y[n][q],tri_z[n][q],c(0),c(1),c(2));
}

void sixdof_geometry::transform(const Eigen::Matrix3d &R, const Eigen::Vector3d &c)
{
	for(n=0; n<tricount; ++n)
	{
        for(int q=0; q<3; q++)
        {
            Eigen::Vector3d point(tri_x0[n][q], tri_y0[n][q], tri_z0[n][q]);
					
            point = R*point;
        
            tri_x[n][q] = point(0) + c(0);
            tri_y[n][q] = point(1) + c(1);
            tri_z[n][q] = point(2) + c(2);
        }
	}
}

void sixdof_geometry::rotation_tri(double phi_, double theta_, double psi_, 
                                   double &xvec, double &yvec, double &zvec, 
                                   const double& x0, const double& y0, const double& z0)
{
	// Distance to origin
    double dx = xvec - x0;
    double dy = yvec - y0;
    double dz = zvec - z0;
	
	// Rotation using Goldstein page 603 (but there is wrong result)
    xvec = dx*(cos(psi_)*cos(theta_)) + dy*(cos(theta_)*sin(psi_)) - dz*sin(theta_);
    yvec = dx*(cos(psi_)*sin(phi_)*sin(theta_)-cos(phi_)*sin(psi_)) + dy*(cos(phi_)*cos(psi_)+sin(phi_)*sin(psi_)*sin(theta_)) + dz*(cos(theta_)*sin(phi_));
    zvec = dx*(sin(phi_)*sin(psi_)+cos(phi_)*cos(psi_)*sin(theta_)) + dy*(cos(phi_)*sin(psi_)*sin(theta_)-cos(psi_)*sin(phi_)) + dz*(cos(phi_)*cos(theta_));
	
	// Moving back
    xvec += x0;
    yvec += y0;
    zvec += z0;
}

// ---------------------------------------------------------------------------------------------
// properties
// ---------------------------------------------------------------------------------------------

double sixdof_geometry::volume() const
{
    double x1, x2, x3, y1, y2, y3, z1, z2, z3;
    double V=0.0;
    
    for (int n = 0; n < tricount; ++n)
    {
        x1 = tri_x[n][0];
        x2 = tri_x[n][1];
        x3 = tri_x[n][2];
        
        y1 = tri_y[n][0];
        y2 = tri_y[n][1];
        y3 = tri_y[n][2];
        
        z1 = tri_z[n][0];
        z2 = tri_z[n][1];
        z3 = tri_z[n][2];  
    
        V += (1.0/6.0)*(-x3*y2*z1 + x2*y3*z1 + x3*y1*z2 - x1*y3*z2 - x2*y1*z3 + x1*y2*z3);
    }
    
    return V;
}
