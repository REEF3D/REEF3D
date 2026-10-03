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

#include"nhflow_amr_6dof.h"
#include"6DOF_nhflow.h"
#include"6DOF_obj.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

nhflow_amr_6dof::nhflow_amr_6dof(lexer *p, sixdof_nhflow *bb, double ds) : b(bb), pp(p), IO(nullptr), CL(nullptr), CR(nullptr), dsm(ds)
{
    n7 = p->imax*p->jmax*(p->kmax+2);

    if(b!=nullptr)
    {
    p->Iarray(IO,n7);
    p->Iarray(CL,n7);
    p->Iarray(CR,n7);
    }
}

nhflow_amr_6dof::~nhflow_amr_6dof()
{
    if(b!=nullptr)
    {
    pp->del_Iarray(IO,n7);
    pp->del_Iarray(CL,n7);
    pp->del_Iarray(CR,n7);
    }
}

// the hull on the grid; the bounds of the ray cast (the subdomain on level 0) are the computed
// cells of the patch
void nhflow_amr_6dof::ray_cast(lexer *p, fdm_nhf *d, ghostcell *pgc, int nb)
{
    const double ox=p->originx, oy=p->originy, ex=p->endx, ey=p->endy;
    p->originx = p->XN[marge];
    p->endx = p->XN[p->knox+marge];
    p->originy = p->YN[marge];
    p->endy = p->YN[p->knoy+marge];

    b->object(nb)->ray_cast_nhflow_grid(p,d,pgc,IO,CL,CR,dsm);

    p->originx=ox; p->originy=oy; p->endx=ex; p->endy=ey;
}

void nhflow_amr_6dof::body(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    if(b==nullptr)
    return;

    for(int nb=0; nb<b->objects(); ++nb)
    ray_cast(p,d,pgc,nb);

    // solid flags as nhflow_forcing::forcing for a floating body (X 16 0)
    LOOP
    p->DF[IJK] = (d->FB[IJK]<0.0) ? -2 : 1;
    pgc->startintV(p,p->DF,1);

    FLOOP
    p->DFF[FIJK] = 1;

    FLOOP
    if(k>0 && k<=p->knoz)
    if(p->DF[IJKm1]<0 && p->DF[IJK]<0)
    p->DFF[FIJK] = -1;

    pgc->startintVF(p,p->DFF,1);
}

// stage iter, before the projection: the body was advanced on level 0
void nhflow_amr_6dof::start_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, int iter,
                                   double *U, double *V, double *W, double *FX, double *FY, double *FZ, slice &WL, slice &fe, bool finalize)
{
    if(b==nullptr)
    return;

    for(int nb=0; nb<b->objects(); ++nb)
    {
        sixdof_obj *o = b->object(nb);
        ray_cast(p,d,pgc,nb);
        o->update_forcing_nhflow(p,d,pgc,d->U,d->V,d->W,FX,FY,FZ,WL,fe,iter);
    }
}

// stage iter, after the projection: the rigid-body velocity again (the loads come from level 0)
void nhflow_amr_6dof::reforce_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, int iter,
                                     double *U, double *V, double *W, double *FX, double *FY, double *FZ, slice &WL, slice &fe, bool finalize)
{
    if(b==nullptr)
    return;

    for(int nb=0; nb<b->objects(); ++nb)
    b->object(nb)->update_forcing_nhflow(p,d,pgc,d->U,d->V,d->W,FX,FY,FZ,WL,fe,iter);
}
