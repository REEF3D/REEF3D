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

//  FNPF mesh refinement in checkpoints (P 48, I 45): the arrays of fnpf_ckpt_fields per patch, written
//  and read by the core (reefamr_ckpt.cpp).

#include"fnpf_amr.h"
#include"reefckpt.h"
#include"fnpf_ckpt_fields.h"

using namespace std;

void fnpf_amr::ckpt_write(reefckpt &ck)
{
    ckpt_write_layout(ck);
    ckpt_write_patches(ck,"fnpf_amr",[&](reefamr_patch *q, const ckpt_visit &v)
    {
        fnpf_ckpt_fields(q->pp,FP(q)->c,[&](const char *n, char k, double *V){ v(n,k,V,nullptr); },
                                         [&](const char *n, char k, int *V){ v(n,k,nullptr,V); });
    });
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

    return ckpt_read_patches(ck,"fnpf_amr",[&](reefamr_patch *q, const ckpt_visit &v)
    {
        fnpf_ckpt_fields(q->pp,FP(q)->c,[&](const char *n, char k, double *V){ v(n,k,V,nullptr); },
                                         [&](const char *n, char k, int *V){ v(n,k,nullptr,V); });
    },msg);
}
