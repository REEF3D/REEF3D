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

#include"sflow_amr.h"
#include"hydro_sflow.h"
#include"sflow_momentum_SI.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"

// REEF3D::HYDRO (R 10 - R 50) with the mesh refinement and on multi-block grids.
//
// Every grid has a hydro_sflow of its own (level 0: fdm2D::hydro of the native grid with the inputs of all grids;
// the patches and the blocks 1..: created here, they read the inputs of level 0).  Per step:
//   - hydro_prepare (before the momentum step): the rain records and discharges at the start of the step (level 0),
//     then S and R of every grid at its own cells (the rasters at the cell centres, the point sources in the cells
//     of the grid); the rain rate of the budget from the leaf cells (rsi_leaf)
//   - the kernels of every grid add S in their continuity (sflow_eta::update, sflow_momentum_SI); on the cells a
//     finer patch covers the restriction replaces the result
//   - hydro_sinks (after the momentum step): the losses on every grid, the volumes and the storage of the budget
//     over the leaf cells; the semi-implicit water column of the patches (water_column, level 0 in sflow_f)
// At a regrid a fresh patch takes the cumulative infiltration and rain of the old patch where it had one, else of
// its parent (constant in the parent cell).  G 7 1: the sources hold over the level-0 step, the losses act once
// after it with its time step.

namespace
{
typedef reefamr_comms_off comms_off;
}

hydro_sflow* sflow_amr::hy(int g)
{
    return (g<0) ? b0->hydro : SP(g)->b->hydro;
}

void sflow_amr::hydro_attach(int id, ghostcell *pgc)
{
    sflow_amr_patch *c = SP(id);

    if(c->b->hydro!=nullptr || b0->hydro==nullptr)
    return;

    c->b->hydro = new hydro_sflow(c->pp,c->b,pgc,b0->hydro);
    c->b->hydro->grid_ini(c->pp,c->b);
}

void sflow_amr::hydro_prepare(lexer *p, fdm2D *b, ghostcell *pgc)
{
    hydro_sflow *h0 = b0->hydro;

    h0->inputs_step(p,pgc);
    h0->grid_prepare(p);

    double q = h0->grid_rain(p,[&](int ii, int jj){ return rsi_leaf(-1,ii,jj); });

    {
    comms_off guard(pgc);
    for(int l=0; l<=maxlev; ++l)
    for(int id : lev[l])
    {
        hydro_attach(id,pgc);
        sflow_amr_patch *c = SP(id);
        c->b->hydro->grid_prepare(c->pp);
        q += c->b->hydro->grid_rain(c->pp,[&](int ii, int jj){ return rsi_leaf(id,ii,jj); });
    }
    }

    h0->step_rain(pgc->globalsum(q));

    // the storage at the start over the leaf cells of all grids (level 0 alone at hydro_sflow::ini, before the patches)
    if(!h0->store0_leaf)
    {
    h0->Vstore0 = h0->Vstore = hydro_storage(p,pgc);
    h0->store0_leaf = true;
    }
}

double sflow_amr::hydro_storage(lexer *p, ghostcell *pgc)
{
    double v = b0->hydro->grid_storage(p,b0,[&](int ii, int jj){ return rsi_leaf(-1,ii,jj); });

    for(int l=0; l<=maxlev; ++l)
    for(int id : lev[l])
    {
        sflow_amr_patch *c = SP(id);
        if(c->b->hydro!=nullptr)
        v += c->b->hydro->grid_storage(c->pp,c->b,[&](int ii, int jj){ return rsi_leaf(id,ii,jj); });
    }

    return pgc->globalsum(v);
}

void sflow_amr::hydro_sinks(lexer *p, fdm2D *b, ghostcell *pgc)
{
    hydro_sflow *h0 = b0->hydro;
    const double dt = p->dt;
    double v[3] = {0.0,0.0,0.0};
    double store = 0.0;

    h0->grid_sinks(p,b0,dt,[&](int ii, int jj){ return rsi_leaf(-1,ii,jj); },v);
    store += h0->grid_storage(p,b0,[&](int ii, int jj){ return rsi_leaf(-1,ii,jj); });

    {
    comms_off guard(pgc);
    for(int l=0; l<=maxlev; ++l)
    for(int id : lev[l])
    {
        sflow_amr_patch *c = SP(id);
        if(c->b->hydro==nullptr)
        continue;
        auto leaf = [&](int ii, int jj){ return rsi_leaf(id,ii,jj); };
        c->b->hydro->grid_sinks(c->pp,c->b,dt,leaf,v);
        store += c->b->hydro->grid_storage(c->pp,c->b,leaf);
        if(c->psim!=nullptr)
        c->psim->water_column(c->pp,c->b);
    }
    }

    h0->step_volumes(p,pgc,v[0],v[1],v[2],store);
}

void sflow_amr::hydro_regrid(ghostcell *pgc, int l, vector<reefamr_patch*> &oldP)
{
    if(b0->hydro==nullptr)
    return;

    {
    comms_off guard(pgc);
    for(int id : lev[l])
    if(P[id]->fresh)
    {
        delete SP(id)->b->hydro;
        SP(id)->b->hydro = nullptr;
        hydro_attach(id,pgc);
    }
    }

    // the cells an old patch held
    old_run(l,2,8460+l,oldP,
            [&](reefamr_patch *q, int io, int jo, double *v)
            {
                hydro_sflow *h = SP(q)->b->hydro;
                v[0] = (h!=nullptr) ? h->F(io,jo) : 0.0;
                v[1] = (h!=nullptr) ? h->Rc(io,jo) : 0.0;
            },
            [&](reefamr_patch *q, int id, int ii, int jj, const double *v)
            {
                hydro_sflow *h = SP(q)->b->hydro;
                h->F(ii,jj) = v[0];
                h->Rc(ii,jj) = v[1];
            });

    // the other cells from the parent, the parent value in all its children
    block_down_if(l,2,8480+l,[](reefamr_patch *q) { return q->fresh; },
                  [&](const reefamr_block &B, int key, double *v)
                  {
                      hydro_sflow *h = hy(B.g);
                      v[0] = (h!=nullptr) ? h->F(B.ic,B.jc) : 0.0;
                      v[1] = (h!=nullptr) ? h->Rc(B.ic,B.jc) : 0.0;
                  },
                  [&](reefamr_patch *q, int id, int k, const double *v)
                  {
                      sflow_amr_patch *c = SP(q);
                      if(!q->fresh || c->blkold[k])
                      return;
                      const int nby = q->ny>>sy;
                      const int i0 = EXT+2*(k/nby), j0 = EY+ychild(k%nby);
                      for(int a=0; a<2; ++a)
                      for(int d=0; d<=sy; ++d)
                      {
                      c->b->hydro->F(i0+a,j0+d) = v[0];
                      c->b->hydro->Rc(i0+a,j0+d) = v[1];
                      }
                  });
}
