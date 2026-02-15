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
Architect: Hans Bihs
--------------------------------------------------------------------*/

#include"6DOF_obj_cfd.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"fieldint.h"
#include"reinidisc.h"

//  The level set fb of a body on the grid of a patch of the CFD mesh refinement (cfd_amr_6dof):
//  as ray_cast and reini_RK2, with the workspace of the patch grid and over the cells around the
//  patch as well (the neighbours of the forcing at the patch edges).

void sixdof_obj_cfd::ray_cast_amr(lexer *p, fdm *a, ghostcell *pgc, fieldint5 &io, fieldint5 &cl, fieldint5 &cr)
{
    const int m = p->margin;      // the cells around the grid (increment::marge is that of the node arrays)
    const int jm = (p->j_dir==1) ? m : 0;

    for(i=-m; i<p->knox+m; ++i)
    for(j=-jm; j<p->knoy+jm; ++j)
    for(k=-m; k<p->knoz+m; ++k)
    {
    io(i,j,k)=1;
    a->fb(i,j,k)=1.0e8;
    }

    // the cells of the field, those outside the domain excepted (geo_raycast_cart: valid)
    geo_cart g = geo_raycast::interior(p);
    g.is = -m;  g.ie = p->knox+m;
    g.js = -jm; g.je = p->knoy+jm;
    g.ks = -m;  g.ke = p->knoz+m;
    g.interior = false;

    for(rayiter=0; rayiter<2; ++rayiter)
    {
        for(int qn=0;qn<entity_sum;++qn)
        {
            if(rayiter==0)
            {
            georay.cart_io(p,0,tri_x,tri_y,tri_z,tstart[qn],tend[qn],g,io.data(),cl.data(),cr.data());

            if(p->j_dir==1)
            georay.cart_io(p,1,tri_x,tri_y,tri_z,tstart[qn],tend[qn],g,io.data(),cl.data(),cr.data());

            georay.cart_io(p,2,tri_x,tri_y,tri_z,tstart[qn],tend[qn],g,io.data(),cl.data(),cr.data());
            }

            if(rayiter==1 && p->X188==1)
            {
            georay.cart_dist(p,0,tri_x,tri_y,tri_z,tstart[qn],tend[qn],g,io.data(),a->fb.data());

            if(p->j_dir==1)
            georay.cart_dist(p,1,tri_x,tri_y,tri_z,tstart[qn],tend[qn],g,io.data(),a->fb.data());

            georay.cart_dist(p,2,tri_x,tri_y,tri_z,tstart[qn],tend[qn],g,io.data(),a->fb.data());
            }

            if(rayiter==1 && p->X188==2)
            georay.cart_vertexdist(p,tri_x,tri_y,tri_z,tstart[qn],tend[qn],a->fb.data());
        }
    }

    const double fb_eps = 1.0e-10*p->DXM;

    for(i=-m; i<p->knox+m; ++i)
    for(j=-jm; j<p->knoy+jm; ++j)
    for(k=-m; k<p->knoz+m; ++k)
    {
        if(io(i,j,k)==-1)
        a->fb(i,j,k)=-fabs(a->fb(i,j,k));

        if(io(i,j,k)==1)
        a->fb(i,j,k)=fabs(a->fb(i,j,k));

        if(fabs(a->fb(i,j,k))<fb_eps)
        a->fb(i,j,k)=0.0;

        a->fb(i,j,k) = MAX(MIN(a->fb(i,j,k),10.0*p->DXM),-10.0*p->DXM);
    }

    // the outermost layer (no neighbour for the distance check): from the next layer inside
    for(i=-m; i<p->knox+m; ++i)
    for(j=-jm; j<p->knoy+jm; ++j)
    for(k=-m; k<p->knoz+m; ++k)
    {
        const int ii = MAX(MIN(i,p->knox+m-2),1-m);
        const int jj = (p->j_dir==1) ? MAX(MIN(j,p->knoy+m-2),1-m) : j;
        const int kk = MAX(MIN(k,p->knoz+m-2),1-m);
        if(ii!=i || jj!=j || kk!=k)
        a->fb(i,j,k) = a->fb(ii,jj,kk);
    }

    pgc->start4a(p,a->fb,50);
}

// reini_RK2 on the cells of the patch; the cells around it keep the distance of the ray cast
void sixdof_obj_cfd::reini_amr(lexer *p, fdm *a, ghostcell *pgc, field &f, reinidisc *prd, field4a &frk, field4a &LL, field4a &dtt)
{
    LOOP
    {
    if(p->j_dir==0)
    dtt.data()[IJK] = p->F43*MIN(p->DXP[IP],p->DZP[KP]);

    if(p->j_dir==1)
	dtt.data()[IJK] = p->F43*MIN3(p->DXP[IP],p->DYP[JP],p->DZP[KP]);
	}

    const int n = p->imax*p->jmax*p->kmax;
    for(int q=0; q<n; ++q)
    frk.data()[q] = f.data()[q];

    for(int q=0;q<5;++q)
	{
		prd->start(p,a,pgc,f,LL,5);

		LOOP
		frk.data()[IJK] = f.data()[IJK] + dtt.data()[IJK]*LL.data()[IJK];

        pgc->start4a(p,frk,50);

		prd->start(p,a,pgc,frk,LL,5);

		LOOP
		f.data()[IJK] = 0.5*f.data()[IJK] + 0.5*frk.data()[IJK] + 0.5*dtt.data()[IJK]*LL.data()[IJK];

        pgc->start4a(p,f,50);
	}
}

// G 7 1 (cfd_amr): stage iter of the body with the step of the level (p->dt): the loads (Xe .. Ne)
// of the last stage of the finest level, the external forces of the stage; the pose, the trimesh
// and the body velocities.  The patches cast the hull at their stages (cfd_amr_6dof)
void sixdof_obj_cfd::amr_stage_cfd(lexer *p, fdm *a, ghostcell *pgc, int iter)
{
    solve_eqmotion_cfd(p,a,pgc,iter,false);
    quat_matrices(p);
    pose_update(p,pgc);
}
