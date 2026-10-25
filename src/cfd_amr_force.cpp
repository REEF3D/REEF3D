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

#include"cfd_amr.h"
#include"force.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"field.h"
#include<algorithm>

//  The forces on the solids (P 81, force) with the mesh refinement.
//
//  force triangulates the zero level of the solid field in its box (P 81) and integrates the
//  pressure and the shear over the pieces.  With the refinement every grid has its own surface:
//  level 0 and every patch that reaches into the box triangulate their field (the patches the
//  surface cast at their resolution, cfd_amr_geom.cpp), and every piece is integrated on the finest
//  grid at its centre (leaf_level: the tile maps of the levels, the same on every rank): level 0
//  leaves out its pieces under a patch, a patch takes those in its interior that no finer patch
//  covers.  A patch that reaches into a force box has a force object without output files (made
//  with the patch, triangulated when first asked) as long as it lives; a regrid makes the new
//  patches theirs.  The sums of the
//  patches of a rank go into the sums of level 0 before the sums over the ranks (force_calc), so
//  the force file and its columns are those of level 0.  The surface file of the box (P 81, the
//  vtp of REEF3D_SOLID) is that of level 0.

namespace
{
struct pscope5
{
    ghostcell *g;
    fdm *back;
    bool oldc, oldl;
    pscope5(ghostcell *gg, fdm *a, fdm *a0) : g(gg), back(a0)
    {
        oldc = g->set_comms(false);
        oldl = g->set_local(true);
        g->fdm_update(a);
    }
    ~pscope5()
    {
        g->fdm_update(back);
        g->set_comms(oldc);
        g->set_local(oldl);
    }
};

// the cell of [0,n) whose nodes hold x (N with the margin offset), -1 outside
int node_cell(const vector<double> &N, int n, int m, double x)
{
    if(x<N[m] || x>=N[n+m])
    return -1;
    int a = 0, b = n-1;
    while(a<b)
    {
        const int c = (a+b+1)/2;
        if(N[c+m]<=x)
        a = c;
        else
        b = c-1;
    }
    return a;
}
}

int cfd_amr::leaf_level(double x, double y, double z) const
{
    const double X[3] = {x,y,z};
    int I[3];
    for(int d=0; d<3; ++d)
    {
        if(gn0[d]==1)
        {
            I[d] = 0;
            continue;
        }
        I[d] = node_cell(gx[d],gn0[d],marge,X[d]);
        if(I[d]<0)
        return 0;
    }

    int l = 0;
    while(l<maxlev && covered(l,I))
    {
        for(int d=0; d<3; ++d)
        {
            const int C = I[d]*rr[d];
            I[d] = (rr[d]==2 && X[d]>=node(l+1,d,C+1)) ? C+1 : C;
        }
        ++l;
    }
    return l;
}

void cfd_amr::solid_forces(int ID, double *F) const
{
    const_cast<cfd_amr*>(this)->patch_forces(ID,F);
}

// the force objects of patch c (patch_objects): one per force box it reaches into.  They are made with
// the patch, not when first asked: their work fields are fields of the patch lexer and so in its
// checkpoint (ckpt_registry, in the order of construction), the same on a restart
void cfd_amr::force_objects(cfd_amr_patch *c)
{
    lexer *p = p0;
    lexer *pp = c->pp;
    const int m = marge;
    if(!solids || p->P81<=0)
    return;

    const double box[6] = {pp->XN[m],pp->XN[pp->knox+m],pp->YN[m],pp->YN[pp->knoy+m],pp->ZN[m],pp->ZN[pp->knoz+m]};
    // in 2D the patch spans the domain in y, as level 0
    double b[6] = {box[0],box[1],box[2],box[3],box[4],box[5]};
    if(pp->j_dir==0)
    {
        b[2] = -1.0e20;
        b[3] = 1.0e20;
    }

    c->forces.assign(p->P81,nullptr);
    c->forces_tri.assign(p->P81,0);
    for(int ID=0; ID<p->P81; ++ID)
    {
        // the interior of the patch against the box of the force
        if(box[1]<p->P81_xs[ID] || box[0]>p->P81_xe[ID] || box[5]<p->P81_zs[ID] || box[4]>p->P81_ze[ID])
        continue;
        if(pp->j_dir==1 && (box[3]<p->P81_ys[ID] || box[2]>p->P81_ye[ID]))
        continue;

        c->forces[ID] = new force(pp,c->a,pgc0,ID,false);
        c->forces[ID]->amr_patch(this,c->lev,b);
    }
}

void cfd_amr::patch_forces(int ID, double *F)
{
    lexer *p = p0;
    if(!solids || ID<0 || ID>=p->P81)
    return;

    for(int id=0; id<npat; ++id)
    {
        cfd_amr_patch *c = CP(id);
        lexer *pp = c->pp;

        if(ID>=(int)c->forces.size() || c->forces[ID]==nullptr)
        continue;

        pscope5 ps(pgc0,c->a,a0);

        // the triangulation of the surface of the patch (its solid field is static)
        if(c->forces_tri[ID]==0)
        {
            // the triangulation reads the solid field of the cells around the patch: those outside the
            // domain (the boundary values of the patch; in 2D the cells across y) take the nearest
            // cell in the domain for it, so that no surface is made at the domain boundary (as the
            // boundary values of level 0)
            field &sf = c->a->solid;
            vector<double> keep(sf.V,sf.V+(size_t)pp->imax*pp->jmax*pp->kmax);
            auto clampd = [](int a, int lo, int hi) { return a<lo ? lo : (a>hi ? hi : a); };
            for(int i=pp->imin; i<pp->imax+pp->imin; ++i)
            for(int j=pp->jmin; j<pp->jmax+pp->jmin; ++j)
            for(int k=pp->kmin; k<pp->kmax+pp->kmin; ++k)
            {
                const int ii = clampd(i,-pp->origin_i,pp->gknox-1-pp->origin_i);
                const int jj = clampd(j,-pp->origin_j,pp->gknoy-1-pp->origin_j);
                const int kk = clampd(k,-pp->origin_k,pp->gknoz-1-pp->origin_k);
                if(ii!=i || jj!=j || kk!=k)
                sf(i,j,k) = keep[(size_t)(ii-pp->imin)*pp->jmax*pp->kmax + (size_t)(jj-pp->jmin)*pp->kmax + (kk-pp->kmin)];
            }

            c->forces[ID]->ini(pp,c->a,pgc0,false);
            c->forces_tri[ID] = 1;

            std::copy(keep.begin(),keep.end(),sf.V);
        }

        c->forces[ID]->patch_sum(pp,c->a,F);
    }
}
