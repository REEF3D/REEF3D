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
	cout<<"Surface triangles: "<<tricount<<endl;
    
    if(p->mpirank==0)
    print_vtp(p,2);
    
}

void nhflow_geometry::objects_allocate_vrans(lexer *p, ghostcell *pgc)
{
    double U,ds,phi,r,snum,trisum;
    
    entity_sum = p->B210 + p->B212 + p->B213 + p->B214 + p->B215 + p->B217 + p->B218 + p->B219;
	tricount=0;
    trisum=0;
    
    
    // ----------------------
    // box
    box_num = p->B210;
    
    p->Darray(box_xs,box_num);
    p->Darray(box_xe,box_num);
    
    p->Darray(box_ys,box_num);
    p->Darray(box_ye,box_num);
    
    p->Darray(box_zs,box_num);
    p->Darray(box_ze,box_num);
    
    for(n=0; n<p->B210;++n)
	{
    trisum+=12;
    
    box_xs[n] = p->B210_xs[n];
    box_xe[n] = p->B210_xe[n];
    
    box_ys[n] = p->B210_ys[n];
    box_ye[n] = p->B210_ye[n];
    
    box_zs[n] = p->B210_zs[n];
    box_ze[n] = p->B210_ze[n];
    }
    
    // ----------------------
    // cylinder_y
    cyly_num = p->B212;
    
    p->Darray(cyly_xc,cyly_num);
    p->Darray(cyly_zc,cyly_num);
    
    p->Darray(cyly_ys,cyly_num);
    p->Darray(cyly_ye,cyly_num);
    
    p->Darray(cyly_r,cyly_num);

    
    for(n=0; n<p->B212;++n)
	{
	r = p->B212_r[n];
	U = 2.0*PI*r;
	ds = 0.5*(DSM);
    snum = MAX(int(U/ds), 40);
	trisum+=6*snum;
    
    cyly_xc[n] = p->B212_xc[n];
    cyly_zc[n] = p->B212_zc[n];
    
    cyly_ys[n] = p->B212_ys[n];
    cyly_ye[n] = p->B212_ye[n];
    
    cyly_r[n] = p->B212_r[n];
	}
    
    // ----------------------
    // cylinder_z
    cylz_num = p->B213;
    
    p->Darray(cylz_xc,cylz_num);
    p->Darray(cylz_yc,cylz_num);
    
    p->Darray(cylz_zs,cylz_num);
    p->Darray(cylz_ze,cylz_num);
    
    p->Darray(cylz_r,cylz_num);
    
    for(n=0; n<p->B213;++n)
	{
	r = p->B213_r[n];
	U = 2.0*PI*r;
	ds = 0.5*(DSM);
    snum = MAX(int(U/ds), 40);
	trisum+=6*snum;
    
    cylz_xc[n] = p->B213_xc[n];
    cylz_yc[n] = p->B213_yc[n];
    
    cylz_zs[n] = p->B213_zs[n];
    cylz_ze[n] = p->B213_ze[n];
    
    cylz_r[n] = p->B213_r[n];
	}
    
    // ----------------------
    // cylinder_member
    jacket_num = p->B214;
    
    p->Darray(jacket_xm1,jacket_num);
    p->Darray(jacket_ym1,jacket_num);
    p->Darray(jacket_zm1,jacket_num);
    p->Darray(jacket_r1,jacket_num);
    
    p->Darray(jacket_xm2,jacket_num);
    p->Darray(jacket_ym2,jacket_num);
    p->Darray(jacket_zm2,jacket_num);
    p->Darray(jacket_r2,jacket_num);
    
    for(n=0; n<p->B214;++n)
	{
	r = MAX(p->B214_r1[n],p->B214_r2[n]);
	U = 2.0*PI*r;
	ds = 0.5*(DSM);
    snum = MAX(int(U/ds), 40);
	trisum+=6*snum;
    
    jacket_xm1[n] = p->B214_xm1[n];
    jacket_ym1[n] = p->B214_ym1[n];
    jacket_zm1[n] = p->B214_zm1[n];
    jacket_r1[n] = p->B214_r1[n];
    
    jacket_xm2[n] = p->B214_xm2[n];
    jacket_ym2[n] = p->B214_ym2[n];
    jacket_zm2[n] = p->B214_zm2[n];
    jacket_r2[n] = p->B214_r2[n];
	}
    
    // ----------------------
    // sphere
    sphere_num = p->B215;
    
    p->Darray(sphere_xm,sphere_num);
    p->Darray(sphere_ym,sphere_num);
    p->Darray(sphere_zm,sphere_num);
    p->Darray(sphere_r,sphere_num);
    
    for(n=0; n<p->B215;++n)
    {
	r = p->B215_r[n];
	U = 2.0*PI*r;
	ds = 0.85*(DSM);
    snum = MAX(int(U/ds), 40);
    trisum+=snum*snum*2;
    
    sphere_xm[n] = p->B215_xm[n];
    sphere_ym[n] = p->B215_ym[n];
    sphere_zm[n] = p->B215_zm[n];
    sphere_r[n] = p->B215_r[n];
    }
    
    // ----------------------
    // wedge
    wedgex_num = p->B217;
    
    p->Darray(wedgex_xs,wedgex_num);
    p->Darray(wedgex_xe,wedgex_num);
    
    p->Darray(wedgex_ys,wedgex_num);
    p->Darray(wedgex_ye,wedgex_num);
    
    p->Darray(wedgex_zs,wedgex_num);
    p->Darray(wedgex_ze,wedgex_num);
    
    for(n=0; n<p->B217;++n)
    {
    trisum+=8;
    
    wedgex_xs[n] = p->B217_xs[n];
    wedgex_xe[n] = p->B217_xe[n];
    
    wedgex_ys[n] = p->B217_ys[n];
    wedgex_ye[n] = p->B217_ye[n];
    
    wedgex_zs[n] = p->B217_zs[n];
    wedgex_ze[n] = p->B217_ze[n];
    }
    
    // ----------------------
    // wedge
    wedgey_num = p->B218;
    
    p->Darray(wedgey_xs,wedgey_num);
    p->Darray(wedgey_xe,wedgey_num);
    
    p->Darray(wedgey_ys,wedgey_num);
    p->Darray(wedgey_ye,wedgey_num);
    
    p->Darray(wedgey_zs,wedgey_num);
    p->Darray(wedgey_ze,wedgey_num);
    
    for(n=0; n<p->B218;++n)
    {
    trisum+=8;
    
    wedgey_xs[n] = p->B218_xs[n];
    wedgey_xe[n] = p->B218_xe[n];
    
    wedgey_ys[n] = p->B218_ys[n];
    wedgey_ye[n] = p->B218_ye[n];
    
    wedgey_zs[n] = p->B218_zs[n];
    wedgey_ze[n] = p->B218_ze[n];
    }
    
    // ----------------------
    // wedge
    wedgez_num = p->B219;
    
    p->Darray(wedgez_xs,wedgez_num);
    p->Darray(wedgez_xe,wedgez_num);
    
    p->Darray(wedgez_ys,wedgez_num);
    p->Darray(wedgez_ye,wedgez_num);
    
    p->Darray(wedgez_zs,wedgez_num);
    p->Darray(wedgez_ze,wedgez_num);
    
    for(n=0; n<p->B219;++n)
    {
    trisum+=8;
    
    wedgez_xs[n] = p->B219_xs[n];
    wedgez_xe[n] = p->B219_xe[n];
    
    wedgez_ys[n] = p->B219_ys[n];
    wedgez_ye[n] = p->B219_ye[n];
    
    wedgez_zs[n] = p->B219_zs[n];
    wedgez_ze[n] = p->B219_ze[n];
    }
    
    // ----------------------
    // STL
    if(p->B230==1)
    entity_sum+=1;

    p->Darray(tri_x,trisum,3);
	p->Darray(tri_y,trisum,3);
	p->Darray(tri_z,trisum,3); 
    
	p->Iarray(tstart,entity_sum);
	p->Iarray(tend,entity_sum);
}
