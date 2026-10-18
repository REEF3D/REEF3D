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

// true if no vertex coordinate of the bodies is further than tol from the triangles of the last cast
// on this grid (tol 0: the same triangles); otherwise the current triangles are kept for the next
// check (tol < 0: always)
bool cfd_amr_6dof::unchanged(double tol)
{
    const int nob = b->objects();
    bool same = ((int)last.size()==nob);

    for(int nb=0; nb<nob && same; ++nb)
    {
        sixdof_obj_cfd *o = b->object(nb);
        const int n = o->amr_tricount();
        const vector<double> &v = last[nb];
        if((int)v.size()!=9*n)
        {
            same = false;
            break;
        }
        size_t c = 0;
        for(int d=0; d<3 && same; ++d)
        {
            double **t = o->amr_tri(d);
            for(int q=0; q<n && same; ++q)
            for(int e=0; e<3; ++e)
            {
                if(!(fabs(v[c]-t[q][e])<=tol))
                same = false;
                ++c;
            }
        }
    }

    if(!same)
    {
        last.assign(nob,vector<double>());
        for(int nb=0; nb<nob; ++nb)
        {
            sixdof_obj_cfd *o = b->object(nb);
            const int n = o->amr_tricount();
            vector<double> &v = last[nb];
            v.resize(9*n);
            size_t c = 0;
            for(int d=0; d<3; ++d)
            {
                double **t = o->amr_tri(d);
                for(int q=0; q<n; ++q)
                for(int e=0; e<3; ++e)
                v[c++] = t[q][e];
            }
        }
    }
    return same;
}

void cfd_amr_6dof::body(lexer *p, fdm *a, ghostcell *pgc)
{
    unchanged(-1.0);

    for(int nb=0; nb<b->objects(); ++nb)
    cast(p,a,pgc,nb);

    pgc->gcdf_update(p,a);
}

// stage iter: the bodies are at the time of the stage (moved on level 0)
void cfd_amr_6dof::start_cfd(lexer *p, fdm *a, ghostcell *pgc, int iter, field &uvel, field &vvel, field &wvel,
                             field &fx, field &fy, field &fz, bool finalize)
{
    reset_fbh(p,a,pgc);

    // the bodies with the triangles of the last cast keep fb (a cast of one body overwrites fb, so
    // all are cast again when one of them has moved); G 42: up to that many cells off the last cast
    const bool keep = unchanged(MAX(p->G42,0.0)*p->DXM);

    for(int nb=0; nb<b->objects(); ++nb)
    {
        if(!keep)
        cast(p,a,pgc,nb);
        b->object(nb)->update_forcing(p,a,pgc,uvel,vvel,wvel,fx,fy,fz,iter);
    }

    pgc->gcdf_update(p,a);
}
