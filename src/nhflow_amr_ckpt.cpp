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

//  NHFLOW mesh refinement in checkpoints (P 48, I 45): the arrays of nhflow_ckpt_fields and of the
//  turbulence model per patch (the blocks 1.. of a multi-block grid as well), written and read by the
//  core (reefamr_ckpt.cpp).

#include"nhflow_amr.h"
#include"reefckpt.h"
#include"nhflow_ckpt_fields.h"
#include"nhflow_turbulence.h"
#include"nhflow_momentum_func.h"

using namespace std;

namespace
{
void patch_fields(nhflow_amr_patch *c, const reefamr::ckpt_visit &v)
{
    nhflow_ckpt_fields(c->pp,c->d,[&](const char *n, char k, double *V){ v(n,k,V,nullptr); },
                                   [&](const char *n, char k, int *V){ v(n,k,nullptr,V); });
    if(c->pmom!=nullptr)
    c->pmom->ckpt_fields([&](const char *n, char k, double *V){ v(n,k,V,nullptr); },
                         [&](const char *n, char k, int *V){ v(n,k,nullptr,V); });
    if(c->pturb!=nullptr)
    c->pturb->ckpt_fields([&](const char *n, char k, double *V){ v(n,k,V,nullptr); },
                          [&](const char *n, char k, int *V){ v(n,k,nullptr,V); });
}
}

void nhflow_amr::ckpt_write(reefckpt &ck)
{
    ckpt_write_layout(ck);
    ckpt_write_patches(ck,"nhflow_amr",[&](reefamr_patch *q, const ckpt_visit &v){ patch_fields(NP(q),v); });
}

bool nhflow_amr::ckpt_read(reefckpt &ck, string &msg)
{
    // the hierarchy of the checkpoint (a regrid to its tile maps, the state of level 0 read before),
    // then the arrays of its patches
    if(!ckpt_read_layout(ck,msg))
    return false;
    regrid(p0,pgc0,false);

    const bool ok = ckpt_read_patches(ck,"nhflow_amr",[&](reefamr_patch *q, const ckpt_visit &v){ patch_fields(NP(q),v); },msg);
    for(auto q : P)
    if((q->lev>=1 || q->blk>0) && NP(q)->pmom!=nullptr)
    NP(q)->pmom->ckpt_restored();
    return ok;
}
