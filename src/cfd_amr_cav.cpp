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
Author: Ahmet Soydan
--------------------------------------------------------------------*/


#include"cfd_amr.h"
#include"cavitation.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<functional>
#include<vector>

//  Cavitation (W 60) with the mesh refinement (G 7 0): one cavitation object per grid; its vapour
//  fraction is in the registry (cfd_amr_fields.cpp): filled around the patches at every stage,
//  restricted at the end of the stages, new patches from their parent.  The density and viscosity
//  of the cells around the patches and of the patch edges of the composite projection follow from
//  it (rmix, vmix with s = alpha); every grid adds its phase-change terms to its own Poisson rows
//  (pjm_corr::amr_prepare -> cavitation::projection), so the composite projection holds them.
//
//  Step (driver_loop_cfd): the coefficients on level 0 (the driver) and on the patches
//  (cav_coefficients), the step of all grids (step_comp), the vapour fraction on level 0 (the
//  driver) and here on the patches level by level, each with the cells around it filled from the
//  coarser level after that level's update; then the restriction onto the coarser grids and the
//  density and viscosity from the restricted fraction.  The fluxes of alpha through the patch edges
//  are not matched (no flux registers): the vapour volume over the leaf cells is not conserved there.

namespace
{
// patch work: no exchange, rank-local reductions, the fdm of the patch in the ghostcell object (as cfd_amr.cpp)
struct pscope
{
    ghostcell *g;
    fdm *back;
    bool oldc, oldl;
    pscope(ghostcell *gg, fdm *a, fdm *a0) : g(gg), back(a0)
    {
        oldc = g->set_comms(false);
        oldl = g->set_local(true);
        g->fdm_update(a);
    }
    ~pscope()
    {
        g->fdm_update(back);
        g->set_comms(oldc);
        g->set_local(oldl);
    }
};
}

field& cfd_amr::cfield(int g)
{
    if(g<0)
    return cavitation::of(p0)->alpha;
    return CP(g)->pcav->alpha;
}

void cfd_amr::cav_coefficients()
{
    if(!cav)
    return;

    // p_v of level 0 (W 69 changes it at every step)
    const double pv = cavitation::of(p0)->vapour_pressure();
    for(int id=0; id<npat; ++id)
    {
        cfd_amr_patch *c = CP(id);
        pscope ps(pgc0,c->a,a0);
        c->pcav->vapour_pressure(pv);
        c->pcav->coefficients(c->pp,c->a,pgc0);
    }
}

void cfd_amr::cav_step()
{
    if(!cav)
    return;

    vector<fspec> fs;
    fs.push_back({CELL,1,[this](int g) -> field& { return cfield(g); },nullptr});

    // the patches level by level, the cells around them from the coarser level updated before
    for(int l=1; l<=maxlev; ++l)
    {
        fill(l,fs);
        for(int id : lev[l])
        {
            cfd_amr_patch *c = CP(id);
            pscope ps(pgc0,c->a,a0);
            c->pcav->start(c->pp,c->a,pgc0);
        }
    }

    // the cells under the patches: the mean of the children, finest level first
    vector<std::function<field&(int)>> cells, none;
    cells.push_back([this](int g) -> field& { return cfield(g); });
    for(int l=maxlev; l>=1; --l)
    restrict_level(l,cells,none);

    // density and viscosity from the restricted fraction: level 0, the patches (a coarser patch
    // got the restriction of a finer one), and the cells around the patches
    pgc0->start4(p0,cfield(-1),1);
    cavitation::of(p0)->properties(p0,a0,pgc0);
    for(int id=0; id<npat; ++id)
    {
        cfd_amr_patch *c = CP(id);
        pscope ps(pgc0,c->a,a0);
        c->pcav->properties(c->pp,c->a,pgc0);
        margin_rovisc(id);
    }
}
