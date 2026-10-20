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

//  FNPF mesh refinement in checkpoints (P 48, I 45): the arrays of fnpf_ckpt_fields per patch and
//  those of the resolved bodies on the patch (fnpf_6DOF), written and read by the core
//  (reefamr_ckpt.cpp).

#include"fnpf_amr.h"
#include"reefckpt.h"
#include"fnpf_ckpt_fields.h"
#include"fnpf_body.h"

using namespace std;

// the arrays of a patch: those of fnpf_ckpt_fields, then those of the body on the patch
void fnpf_amr::patch_fields(reefamr_patch *q, const ckpt_visit &v)
{
    auto fd = [&](const char *n, char k, double *V){ v(n,k,V,nullptr); };
    auto fi = [&](const char *n, char k, int *V){ v(n,k,nullptr,V); };
    fnpf_ckpt_fields(q->pp,FP(q)->c,fd,fi);
    if(body!=nullptr)
    {
        body->amr_ckpt_fields(q->pp,fd,fi);
        
        // the stage values and tendencies of the patch: the footprint of the body on a patch
        // (amr_surface) reads the cells around its piece before the fills of the stage set them
        fnpf_amr_patch *c = FP(q);
        auto sl = [&](const char *n, slice4 *a){ fd(n,'S',a!=nullptr ? a->V : nullptr); };
        sl("amr_erk1",c->erk1); sl("amr_erk2",c->erk2); sl("amr_frk1",c->frk1); sl("amr_frk2",c->frk2);
        sl("amr_ek",c->ek); sl("amr_fk",c->fk);
    }
}

void fnpf_amr::ckpt_write(reefckpt &ck)
{
    ckpt_write_layout(ck);
    ckpt_write_patches(ck,"fnpf_amr",[&](reefamr_patch *q, const ckpt_visit &v){ patch_fields(q,v); });
}

bool fnpf_amr::ckpt_read(reefckpt &ck, string &msg)
{
    // the hierarchy of the checkpoint (a regrid to its tile maps, the state of level 0 read before),
    // then the arrays of its patches
    if(!ckpt_read_layout(ck,msg))
    return false;
    Se0 = &c0->eta;
    Sf0 = &c0->Fifsf;
    stg = 2;
    regrid(p0,pgc0,false);

    // the body grids of these patches (resolved bodies), their arrays come with the patch
    if(body!=nullptr)
    body->amr_ckpt_grids(p0,pgc0);

    return ckpt_read_patches(ck,"fnpf_amr",[&](reefamr_patch *q, const ckpt_visit &v){ patch_fields(q,v); },msg);
}
