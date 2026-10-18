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
#include"reinidisc_f.h"

//  The level set fb of a body on the grid of a patch of the CFD mesh refinement (cfd_amr_6dof):
//  as ray_cast and reini_RK2, with the workspace of the patch grid and over the cells around the
//  patch as well (the neighbours of the forcing at the patch edges).
//
//  Both work in the box of the body on the grid: the cells within amr_fb_band cells of the bounding
//  box of its triangles.  Outside the box fb is the clamped distance 10 dx of the cells far from the
//  body (the forcing reads fb within 1.1 dx of the body only, Hsolidface).  A small body on a large
//  patch (a propeller on its refined box, the coarser patches around it) costs its own size, not
//  that of the patch.

namespace
{
const int amr_fb_band = 6;      // cells between the bounding box of the triangles and the edge of the box

// [s,e) of the cells of [lo,hi) whose centre lies in [x0,x1]
void amr_fb_range(const double *XP, int lo, int hi, double x0, double x1, int &s, int &e)
{
    s = hi;
    e = lo;
    for(int n=lo; n<hi; ++n)
    if(XP[n+increment::marge]>=x0 && XP[n+increment::marge]<=x1)
    {
        s = MIN(s,n);
        e = MAX(e,n+1);
    }
    if(e<=s)
    s = e = lo;
}
}

void sixdof_obj_cfd::ray_cast_amr(lexer *p, fdm *a, ghostcell *pgc, fieldint5 &io, fieldint5 &cl, fieldint5 &cr, int *box)
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

    // the box of the body: the cells of the field within amr_fb_band cells of its triangles
    double xs[3] = {1.0e20,1.0e20,1.0e20}, xe[3] = {-1.0e20,-1.0e20,-1.0e20};
    for(int qn=0; qn<entity_sum; ++qn)
    for(int n=tstart[qn]; n<tend[qn]; ++n)
    for(int q=0; q<3; ++q)
    {
        xs[0] = MIN(xs[0],tri_x[n][q]);  xe[0] = MAX(xe[0],tri_x[n][q]);
        xs[1] = MIN(xs[1],tri_y[n][q]);  xe[1] = MAX(xe[1],tri_y[n][q]);
        xs[2] = MIN(xs[2],tri_z[n][q]);  xe[2] = MAX(xe[2],tri_z[n][q]);
    }
    const double w = double(amr_fb_band)*p->DXM;
    amr_fb_range(p->XP,-m,p->knox+m,xs[0]-w,xe[0]+w,box[0],box[1]);
    if(p->j_dir==1)
    amr_fb_range(p->YP,-jm,p->knoy+jm,xs[1]-w,xe[1]+w,box[2],box[3]);
    else
    {
    box[2] = 0;
    box[3] = p->knoy;
    }
    amr_fb_range(p->ZP,-m,p->knoz+m,xs[2]-w,xe[2]+w,box[4],box[5]);
    const bool empty = (box[1]<=box[0] || box[3]<=box[2] || box[5]<=box[4]);

    // the cells of the box, those outside the domain excepted (geo_raycast_cart: valid)
    geo_cart g = geo_raycast::interior(p);
    g.is = box[0]; g.ie = box[1];
    g.js = box[2]; g.je = box[3];
    g.ks = box[4]; g.ke = box[5];
    g.interior = false;

    for(rayiter=0; rayiter<2 && !empty; ++rayiter)
    {
        for(int qn=0;qn<entity_sum;++qn)
        {
            if(rayiter==0)
            {
            georay.cart_io(p,0,tri_x,tri_y,tri_z,tstart[qn],tend[qn],g,io.V,cl.V,cr.V);

            if(p->j_dir==1)
            georay.cart_io(p,1,tri_x,tri_y,tri_z,tstart[qn],tend[qn],g,io.V,cl.V,cr.V);

            georay.cart_io(p,2,tri_x,tri_y,tri_z,tstart[qn],tend[qn],g,io.V,cl.V,cr.V);
            }

            if(rayiter==1 && p->X188==1)
            {
            georay.cart_dist(p,0,tri_x,tri_y,tri_z,tstart[qn],tend[qn],g,io.V,a->fb.V);

            if(p->j_dir==1)
            georay.cart_dist(p,1,tri_x,tri_y,tri_z,tstart[qn],tend[qn],g,io.V,a->fb.V);

            georay.cart_dist(p,2,tri_x,tri_y,tri_z,tstart[qn],tend[qn],g,io.V,a->fb.V);
            }

            if(rayiter==1 && p->X188==2)
            georay.cart_vertexdist(p,tri_x,tri_y,tri_z,tstart[qn],tend[qn],a->fb.V);
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

// reini_RK2 on the cells of the patch in the box of the body; the cells around the patch keep the
// distance of the ray cast, the cells outside the box the clamped distance
void sixdof_obj_cfd::reini_amr(lexer *p, fdm *a, ghostcell *pgc, field &f, reinidisc_f *prd, field4a &frk, field4a &LL, field4a &dtt, const int *box)
{
    // the cells of the patch in the box
    const int b[6] = {MAX(box[0],0),MIN(box[1],p->knox),MAX(box[2],0),MIN(box[3],p->knoy),MAX(box[4],0),MIN(box[5],p->knoz)};
    if(b[1]<=b[0] || b[3]<=b[2] || b[5]<=b[4])
    return;

    #define BOXLOOP for(i=b[0]; i<b[1]; ++i) for(j=b[2]; j<b[3]; ++j) for(k=b[4]; k<b[5]; ++k) if(p->flag4[IJK]>0)

    BOXLOOP
    {
    if(p->j_dir==0)
    dtt.V[IJK] = p->F43*MIN(p->DXP[IP],p->DZP[KP]);

    if(p->j_dir==1)
	dtt.V[IJK] = p->F43*MIN3(p->DXP[IP],p->DYP[JP],p->DZP[KP]);
	}

    const int n = p->imax*p->jmax*p->kmax;
    for(int q=0; q<n; ++q)
    frk.V[q] = f.V[q];

    for(int q=0;q<5;++q)
	{
		prd->start_box(p,a,pgc,f,LL,b);

		BOXLOOP
		frk.V[IJK] = f.V[IJK] + dtt.V[IJK]*LL.V[IJK];

        pgc->start4a(p,frk,50);

		prd->start_box(p,a,pgc,frk,LL,b);

		BOXLOOP
		f.V[IJK] = 0.5*f.V[IJK] + 0.5*frk.V[IJK] + 0.5*dtt.V[IJK]*LL.V[IJK];

        pgc->start4a(p,f,50);
	}

    #undef BOXLOOP
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
