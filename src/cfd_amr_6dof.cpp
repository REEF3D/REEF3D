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

#include"cfd_amr_6dof.h"
#include"6DOF_cfd.h"
#include"6DOF_obj_cfd.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"reinidisc_f.h"

cfd_amr_6dof::cfd_amr_6dof(lexer *p, sixdof_cfd *bb) : b(bb), cutl(p), cutr(p), fbio(p), frk1(p), L(p), dt(p)
{
    prdisc = new reinidisc_f(p,true);
}

cfd_amr_6dof::~cfd_amr_6dof()
{
    delete prdisc;
}

// body nb on the grid: ray cast, reinitialisation
void cfd_amr_6dof::cast(lexer *p, fdm *a, ghostcell *pgc, int nb)
{
    sixdof_obj_cfd *o = b->object(nb);
    o->ray_cast_amr(p,a,pgc,fbio,cutl,cutr,box);
    o->reini_amr(p,a,pgc,a->fb,prdisc,frk1,L,dt,box);
}

// the solid Heaviside functions of the stage (sixdof_cfd::setup)
void cfd_amr_6dof::reset_fbh(lexer *p, fdm *a, ghostcell *pgc)
{
    ULOOP
    a->fbh1(i,j,k) = 0.0;

    VLOOP
    a->fbh2(i,j,k) = 0.0;

    WLOOP
    a->fbh3(i,j,k) = 0.0;

    LOOP
    a->fbh4(i,j,k) = 0.0;

    pgc->start1(p,a->fbh1,10);
    pgc->start2(p,a->fbh2,11);
    pgc->start3(p,a->fbh3,12);
    pgc->start4(p,a->fbh4,40);
}

void cfd_amr_6dof::body(lexer *p, fdm *a, ghostcell *pgc)
{
    for(int nb=0; nb<b->objects(); ++nb)
    cast(p,a,pgc,nb);

    pgc->gcdf_update(p,a);
}

// stage iter: the bodies are at the time of the stage (moved on level 0)
void cfd_amr_6dof::start_cfd(lexer *p, fdm *a, ghostcell *pgc, int iter, field &uvel, field &vvel, field &wvel,
                             field &fx, field &fy, field &fz, bool finalize)
{
    reset_fbh(p,a,pgc);

    for(int nb=0; nb<b->objects(); ++nb)
    {
        cast(p,a,pgc,nb);
        b->object(nb)->update_forcing(p,a,pgc,uvel,vvel,wvel,fx,fy,fz,iter);
    }

    pgc->gcdf_update(p,a);
}
