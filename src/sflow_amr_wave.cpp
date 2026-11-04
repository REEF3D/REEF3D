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

#include "fdm2D.h"
#include "ghostcell.h"
#include "lexer.h"
#include "met_forcing.h"
#include "seastate_sflow.h"
#include "sflow_amr.h"
#include "slice4.h"

// REEF3D::SEASTATE coupling (A 750 1 in-process, A 750 2 remote waves) with the mesh refinement.
//
// The wave model is coupled with level 0: it takes the environment from the level-0 fields (the covered cells hold
// the restricted state of the patches) and gives the wave force (Fw, Gw), div(M) of the Stokes transport (A 751 2)
// and Hs, Tp (wave-dependent drag W 207 7) on the level-0 cells.  Every patch has a seastate_sflow of its own grid
// that holds these fields only (no wave model, no long-wave boundary: the patches stay off the domain edges), so the
// kernels of the patch apply them as on level 0 (u_source, v_source in the stages, mass_source in the continuity,
// the ramp A 752 at the time of the patch).
//
// The fields of the patches follow from their parents, coarse to fine, after every change on level 0 (a wave step,
// a receive from the remote waves: seastate_coupling::updates) and after every regrid:
//   - the cells of the patch box: linear in the parent cell with minmod slopes over the fluid neighbours (the mean
//     of the children is the parent value: the force on the parent area is kept; next to dry cells, where the
//     fields are zero, the slopes vanish)
//   - the cells around the box: from the patch of the level that holds them, else linear from the parent
// The level-0 fields first take the ghost cells of the partition (the slopes next to a rank boundary as on one
// rank).  W 207 7: the met forcing of every patch takes its Hs, Tp (met_forcing::waves) and updates its roughness.

namespace {
typedef reefamr_comms_off comms_off;

const int NWF = 5;  // Fw, Gw, divM, Hsw, Tpw

slice4 *wfield(seastate_sflow *w, int n)
{
    if(n == 0)
        return &w->force_x();
    if(n == 1)
        return &w->force_y();
    if(n == 2)
        return &w->div_m();
    if(n == 3)
        return &w->Hsw;
    return &w->Tpw;
}
}  // namespace

seastate_sflow *sflow_amr::wv(int g)
{
    return (g < 0) ? b0->wave : SP(g)->b->wave;
}

void sflow_amr::wave_sync(ghostcell *pgc)
{
    if(b0->wave == nullptr || maxlev < 1)
        return;

    const int u = b0->wave->updates();

    if(u == wave_seen && !wave_dirty)
        return;

    wave_seen = u;
    wave_dirty = false;

    wave_fields(pgc);
}

void sflow_amr::wave_fields(ghostcell *pgc)
{
    // the ghost cells of the partition of the level-0 fields (collective)
    for(int n = 0; n < NWF; ++n)
        pgc->gcsl_start4(p0, *wfield(b0->wave, n), 24);

    // the objects of the patches
    {
        comms_off guard(pgc);
        for(auto q : P)
        {
            sflow_amr_patch *c = SP(q);
            if(c->b->wave == nullptr)
            {
                c->b->wave = new seastate_sflow(c->pp, c->b, pgc);
                if(p0->W207 == 7 && p0->W200 > 0 && c->b->met != nullptr)
                    c->b->met->waves(c->pp, &c->b->wave->Hsw, &c->b->wave->Tpw, &c->b->WL);
            }
        }
    }

    if(patches_total == 0)
        return;

    for(int l = 1; l <= maxlev; ++l)
    {
        // the box: linear in the parent (2 or 4 children per parent cell)
        block_down(l, 4 * NWF, 8280 + l, [&](const reefamr_block &B, int key, double *v)
                   {
                       gh G = grid(B.g);
                       seastate_sflow *w = wv(B.g);
                       for(int n=0; n<NWF; ++n)
                       for(int a=0; a<2; ++a)
                       for(int d=0; d<2; ++d)
                       v[4*n+2*a+d] = (d>sy) ? 0.0 : sed_plin(*wfield(w,n),G.q,B.ic,B.jc,a==0?-1:1,(sy==0) ? 0 : (d==0?-1:1)); }, [&](reefamr_patch *q, int id, int k, const double *v)
                   {
                       const int nby = q->ny>>sy;
                       const int i0 = EXT+2*(k/nby), j0 = EY+ychild(k%nby);
                       seastate_sflow *w = SP(q)->b->wave;
                       for(int n=0; n<NWF; ++n)
                       for(int a=0; a<2; ++a)
                       for(int d=0; d<=sy; ++d)
                       (*wfield(w,n))(i0+a,j0+d) = v[4*n+2*a+d]; });

        // the cells around the box
        fill_run(l, NWF, 8300 + l, [&](const reefamr_fill &f, double *v)
                 {
                     for(int n=0; n<NWF; ++n)
                     {
                         if(f.kind==0)
                         v[n] = (*wfield(wv(f.g),n))(f.si,f.sj);
                         else if(f.kind==1)
                         v[n] = sed_plin(*wfield(wv(f.g),n),grid(f.g).q,f.si,f.sj,f.ox,f.oy);
                         else
                         v[n] = 0.0;
                     } }, [&](reefamr_patch *q, int id, const reefamr_fill &f, const double *w)
                 {
                     seastate_sflow *s = SP(q)->b->wave;
                     for(int n=0; n<NWF; ++n)
                     (*wfield(s,n))(f.di,f.dj) = w[n]; });
    }

    // W 207 7: the roughness of the met forcing of the patches from their sea state
    if(p0->W207 == 7 && p0->W200 > 0)
    {
        comms_off guard(pgc);
        for(auto q : P)
        {
            sflow_amr_patch *c = SP(q);
            if(c->b->met != nullptr)
                c->b->met->waves_update(c->pp);
        }
    }
}
