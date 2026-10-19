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

//  SFLOW mesh refinement in checkpoints (P 48, I 45): the arrays of sflow_ckpt_fields and of the
//  momentum, free-surface, pressure and diffusion objects of every patch (reefamr_ckpt.cpp), the
//  hierarchy of the checkpoint restored by a regrid to its tile maps.

#include"sflow_amr.h"
#include"reefckpt.h"
#include"sflow_ckpt_fields.h"
#include"sflow_momentum_RK3.h"
#include"sflow_fsf.h"
#include"sflow_pressure.h"
#include"sflow_diffusion.h"

using namespace std;

namespace
{
void patch_fields(sflow_amr_patch *c, const reefamr::ckpt_visit &v)
{
    auto fd = [&](const char *n, char k, double *V){ v(n,k,V,nullptr); };
    auto fi = [&](const char *n, char k, int *V){ v(n,k,nullptr,V); };
    sflow_ckpt_fields(c->pp,c->b,fd,fi);
    // shoreline: the restricted surface offsets of the covered cells (used by the stages of the next
    // step) and the parent surfaces of the ghost cells
    if(c->sfc.size()==(size_t)c->pp->imax*c->pp->jmax)
    fd("amr_sfc",'S',c->sfc.data());
    if(c->gsf.size()==(size_t)c->pp->imax*c->pp->jmax)
    fd("amr_gsf",'S',c->gsf.data());
    if(c->pmom!=nullptr)
    c->pmom->ckpt_fields(fd,fi);
    if(c->pfsf!=nullptr)
    c->pfsf->ckpt_fields(fd,fi);
    if(c->ppress!=nullptr)
    c->ppress->ckpt_fields(fd,fi);
    if(c->pdiff!=nullptr)
    c->pdiff->ckpt_fields(fd,fi);
}
}

void sflow_amr::ckpt_write(reefckpt &ck)
{
    ckpt_write_layout(ck);
    ckpt_write_patches(ck,"sflow_amr",[&](reefamr_patch *q, const ckpt_visit &v){ patch_fields(SP(q),v); });
}

bool sflow_amr::ckpt_read(reefckpt &ck, string &msg)
{
    // the hierarchy of the checkpoint (a regrid to its tile maps, the state of level 0 read before),
    // then the arrays of its patches
    if(!ckpt_read_layout(ck,msg))
    return false;
    regrid(p0,pgc0,false);

    return ckpt_read_patches(ck,"sflow_amr",[&](reefamr_patch *q, const ckpt_visit &v){ patch_fields(SP(q),v); },msg);
}
