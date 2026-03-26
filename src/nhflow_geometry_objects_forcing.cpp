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

void nhflow_geometry::objects_create_forcing(lexer *p, ghostcell *pgc)
{
    int qn;

    objects_allocate_forcing(p,pgc);
	
    entity_count=0;
	
	for(qn=0;qn<box_num;++qn)
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
    print_vtp(p,1);

    if(p->mpirank==0)
	cout<<"Surface triangles: "<<tricount<<endl;
}

void nhflow_geometry::objects_allocate_forcing(lexer *p, ghostcell *pgc)
{
    double U,ds,phi,r,snum,trisum;
    
    entity_sum = p->A581 + p->A583 + p->A584 + p->A585 + p->A586 + p->A587 + p->A588 + p->A589;
	tricount=0;
    trisum=0;
    
    
    // ----------------------
    // box
    box_num = p->A581;
    
    p->Darray(box_xs,box_num);
    p->Darray(box_xe,box_num);
    
    p->Darray(box_ys,box_num);
    p->Darray(box_ye,box_num);
    
    p->Darray(box_zs,box_num);
    p->Darray(box_ze,box_num);
    
    for(n=0; n<p->A581;++n)
	{
    trisum+=12;
    
    box_xs[n] = p->A581_xs[n];
    box_xe[n] = p->A581_xe[n];
    
    box_ys[n] = p->A581_ys[n];
    box_ye[n] = p->A581_ye[n];
    
    box_zs[n] = p->A581_zs[n];
    box_ze[n] = p->A581_ze[n];
    }
    
    // ----------------------
    // cylinder_y
    cyly_num = p->A583;
    
    p->Darray(cyly_xc,cyly_num);
    p->Darray(cyly_zc,cyly_num);
    
    p->Darray(cyly_ys,cyly_num);
    p->Darray(cyly_ye,cyly_num);
    
    p->Darray(cyly_r,cyly_num);

    
    for(n=0; n<p->A583;++n)
	{
	r = p->A583_r[n];
	U = 2.0*PI*r;
	ds = 0.5*(DSM);
    snum = MAX(int(U/ds), 40);
	trisum+=6*snum;
    
    cyly_xc[n] = p->A583_xc[n];
    cyly_zc[n] = p->A583_zc[n];
    
    cyly_ys[n] = p->A583_ys[n];
    cyly_ye[n] = p->A583_ye[n];
    
    cyly_r[n] = p->A583_r[n];
	}
    
    // ----------------------
    // cylinder_z
    cylz_num = p->A584;
    
    p->Darray(cylz_xc,cylz_num);
    p->Darray(cylz_yc,cylz_num);
    
    p->Darray(cylz_zs,cylz_num);
    p->Darray(cylz_ze,cylz_num);
    
    p->Darray(cylz_r,cylz_num);
    
    for(n=0; n<p->A584;++n)
	{
	r = p->A584_r[n];
	U = 2.0*PI*r;
	ds = 0.5*(DSM);
    snum = MAX(int(U/ds), 40);
	trisum+=6*snum;
    
    cylz_xc[n] = p->A584_xc[n];
    cylz_yc[n] = p->A584_yc[n];
    
    cylz_zs[n] = p->A584_zs[n];
    cylz_ze[n] = p->A584_ze[n];
    
    cylz_r[n] = p->A584_r[n];
	}
    
    // ----------------------
    // cylinder_member
    jacket_num = p->A585;
    
    p->Darray(jacket_xm1,jacket_num);
    p->Darray(jacket_ym1,jacket_num);
    p->Darray(jacket_zm1,jacket_num);
    p->Darray(jacket_r1,jacket_num);
    
    p->Darray(jacket_xm2,jacket_num);
    p->Darray(jacket_ym2,jacket_num);
    p->Darray(jacket_zm2,jacket_num);
    p->Darray(jacket_r2,jacket_num);
    
    for(n=0; n<p->A585;++n)
	{
	r = MAX(p->A585_r1[n],p->A585_r2[n]);
	U = 2.0*PI*r;
	ds = 0.5*(DSM);
    snum = MAX(int(U/ds), 40);
	trisum+=6*snum;
    
    jacket_xm1[n] = p->A585_xm1[n];
    jacket_ym1[n] = p->A585_ym1[n];
    jacket_zm1[n] = p->A585_zm1[n];
    jacket_r1[n] = p->A585_r1[n];
    
    jacket_xm2[n] = p->A585_xm2[n];
    jacket_ym2[n] = p->A585_ym2[n];
    jacket_zm2[n] = p->A585_zm2[n];
    jacket_r2[n] = p->A585_r2[n];
	}
    
    // ----------------------
    // sphere
    sphere_num = p->A586;
    
    p->Darray(sphere_xm,sphere_num);
    p->Darray(sphere_ym,sphere_num);
    p->Darray(sphere_zm,sphere_num);
    p->Darray(sphere_r,sphere_num);
    
    for(n=0; n<p->A586;++n)
    {
	r = p->A586_r[n];
	U = 2.0*PI*r;
	ds = 0.85*(DSM);
    snum = MAX(int(U/ds), 40);
    trisum+=snum*snum*2;
    
    sphere_xm[n] = p->A586_xm[n];
    sphere_ym[n] = p->A586_ym[n];
    sphere_zm[n] = p->A586_zm[n];
    sphere_r[n] = p->A586_r[n];
    }
    
    // ----------------------
    // wedge
    wedgex_num = p->A587;
    
    p->Darray(wedgex_xs,wedgex_num);
    p->Darray(wedgex_xe,wedgex_num);
    
    p->Darray(wedgex_ys,wedgex_num);
    p->Darray(wedgex_ye,wedgex_num);
    
    p->Darray(wedgex_zs,wedgex_num);
    p->Darray(wedgex_ze,wedgex_num);
    
    for(n=0; n<p->A587;++n)
    {
    trisum+=8;
    
    wedgex_xs[n] = p->A587_xs[n];
    wedgex_xe[n] = p->A587_xe[n];
    
    wedgex_ys[n] = p->A587_ys[n];
    wedgex_ye[n] = p->A587_ye[n];
    
    wedgex_zs[n] = p->A587_zs[n];
    wedgex_ze[n] = p->A587_ze[n];
    }
    
    // ----------------------
    // wedge
    wedgey_num = p->A588;
    
    p->Darray(wedgey_xs,wedgey_num);
    p->Darray(wedgey_xe,wedgey_num);
    
    p->Darray(wedgey_ys,wedgey_num);
    p->Darray(wedgey_ye,wedgey_num);
    
    p->Darray(wedgey_zs,wedgey_num);
    p->Darray(wedgey_ze,wedgey_num);
    
    for(n=0; n<p->A588;++n)
    {
    trisum+=8;
    
    wedgey_xs[n] = p->A588_xs[n];
    wedgey_xe[n] = p->A588_xe[n];
    
    wedgey_ys[n] = p->A588_ys[n];
    wedgey_ye[n] = p->A588_ye[n];
    
    wedgey_zs[n] = p->A588_zs[n];
    wedgey_ze[n] = p->A588_ze[n];
    }
    
    // ----------------------
    // wedge
    wedgez_num = p->A589;
    
    p->Darray(wedgez_xs,wedgez_num);
    p->Darray(wedgez_xe,wedgez_num);
    
    p->Darray(wedgez_ys,wedgez_num);
    p->Darray(wedgez_ye,wedgez_num);
    
    p->Darray(wedgez_zs,wedgez_num);
    p->Darray(wedgez_ze,wedgez_num);
    
    for(n=0; n<p->A589;++n)
    {
    trisum+=8;
    
    wedgez_xs[n] = p->A589_xs[n];
    wedgez_xe[n] = p->A589_xe[n];
    
    wedgez_ys[n] = p->A589_ys[n];
    wedgez_ye[n] = p->A589_ye[n];
    
    wedgez_zs[n] = p->A589_zs[n];
    wedgez_ze[n] = p->A589_ze[n];
    }
    
    // ----------------------
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
