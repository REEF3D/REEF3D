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
//
//  Sediment (S 10, nhflow_amr_sed.cpp): every patch has its own sediment_f.  Its slices (bed, bed load,
//  shear stress, the arrays of the Exner step) are in the registry of the patch lexer (ckpt_registry.h),
//  as those of level 0 under the tag "reg" (driver_ckpt.cpp): the slices and fields of the registry
//  follow the arrays above in the list of the patch.  Its scalars (sediment_f::ckpt_state: the volume,
//  the Exner scalars, the time filter of the bed shear) per patch in the rank blob nhflow_amr_sed with
//  the level and box of the patch; a patch of the restart takes those of the checkpoint patch with
//  the same level and box, else of one of its level (another decomposition), else it keeps its own.

#include"nhflow_amr.h"
#include"reefckpt.h"
#include"nhflow_longrun.h"
#include"nhflow_ckpt_fields.h"
#include"nhflow_turbulence.h"
#include"nhflow_momentum_func.h"
#include"ckpt_registry.h"
#include"sediment_f.h"

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

    // S 10: the slices and fields of the registry of the patch (the sediment objects; the arrays of the
    // flow among them once more), in the order of their construction; raw arrays are not grid arrays
    if(c->sed!=nullptr)
    ckpt_registry::visit(c->pp,[&](const ckpt_registry::entry &e)
    {
        if(e.kind!='S' && e.kind!='C')
        return;
        if(e.isint)
        v("reg",e.kind,nullptr,(int*)e.V);
        else
        v("reg",e.kind,(double*)e.V,nullptr);
    });
}

// S 10: the scalars of the sediment of every patch of this rank (lev, I0, I1, J0, J1, n, state)
void sed_write(const vector<reefamr_patch*> &P, reefckpt &ck)
{
    vector<double> b;
    for(auto q : P)
    {
        nhflow_amr_patch *c = static_cast<nhflow_amr_patch*>(q);
        if((q->lev<1 && q->blk==0) || c->sed==nullptr)
        continue;
        vector<double> sv;
        c->sed->ckpt_state(c->pp,sv);
        const double hd[6] = {double(q->lev),double(q->I0),double(q->I1),double(q->J0),double(q->J1),double(sv.size())};
        b.insert(b.end(),hd,hd+6);
        b.insert(b.end(),sv.begin(),sv.end());
    }
    ck.rank_blob("nhflow_amr_sed",b);
}

bool sed_read(const vector<reefamr_patch*> &P, reefckpt &ck, ghostcell *pgc, string &msg)
{
    struct rec { int lev, I0, I1, J0, J1; vector<double> sv; };
    vector<rec> all;
    for(const auto &bl : ck.get_rank_blobs("nhflow_amr_sed"))
    {
        size_t q=0;
        while(q+6<=bl.size())
        {
            rec r{int(bl[q]),int(bl[q+1]),int(bl[q+2]),int(bl[q+3]),int(bl[q+4]),{}};
            const size_t n = size_t(bl[q+5]);
            q += 6;
            if(q+n>bl.size())
            {
                msg = "the sediment of the mesh refinement in the checkpoint is incomplete";
                return false;
            }
            r.sv.assign(bl.begin()+q,bl.begin()+q+n);
            q += n;
            all.push_back(std::move(r));
        }
    }

    for(auto q : P)
    {
        nhflow_amr_patch *c = static_cast<nhflow_amr_patch*>(q);
        if((q->lev<1 && q->blk==0) || c->sed==nullptr)
        continue;
        const rec *src = nullptr;
        for(const auto &r : all)
        if(r.lev==q->lev && r.I0==q->I0 && r.I1==q->I1 && r.J0==q->J0 && r.J1==q->J1)
        {
            src = &r;
            break;
        }
        if(src==nullptr)
        for(const auto &r : all)
        if(r.lev==q->lev)
        {
            src = &r;
            break;
        }
        if(src!=nullptr && !c->sed->ckpt_restore(c->pp,pgc,src->sv))
        {
            msg = "the sediment state of a patch in the checkpoint does not fit this run (S 10)";
            return false;
        }
    }
    return true;
}
}

void nhflow_amr::ckpt_write(reefckpt &ck)
{
    ckpt_write_layout(ck);
    ckpt_write_patches(ck,"nhflow_amr",[&](reefamr_patch *q, const ckpt_visit &v){ patch_fields(NP(q),v); });
    // the number of the layout (the log of the refinement, REEF3D_NHFLOW_AMR_log.dat)
    ck.blob("nhflow_amr_layout_id",vector<double>(1,double(layout_id)));

    if(p0->S10>0)
    sed_write(P,ck);
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

    if(ok && p0->S10>0 && !sed_read(P,ck,pgc0,msg))
    return false;

    // the layout continues the numbering of the run (a checkpoint without it: counted from the restart)
    if(ck.has("nhflow_amr_layout_id"))
    {
        const vector<double> b = ck.get_blob("nhflow_amr_layout_id");
        if(!b.empty())
        layout_id = int(b[0]);
    }
    return ok;
}

// long-run statistics (P 160): the sums of the blocks with the patches (tag "nhflow_longrun"); the patches
// that are no grid of the record (refined patches) carry no arrays
void nhflow_amr::longrun_ckpt_write(reefckpt &ck, nhflow_longrun *lr)
{
    ckpt_write_patches(ck,"nhflow_longrun",[&](reefamr_patch *q, const ckpt_visit &v)
    {
        lr->ckpt_fields(NP(q)->pp,[&](const char *n, char k, double *V){ v(n,k,V,nullptr); });
    });
}

bool nhflow_amr::longrun_ckpt_read(reefckpt &ck, nhflow_longrun *lr, string &msg)
{
    return ckpt_read_patches(ck,"nhflow_longrun",[&](reefamr_patch *q, const ckpt_visit &v)
    {
        lr->ckpt_fields(NP(q)->pp,[&](const char *n, char k, double *V){ v(n,k,V,nullptr); });
    },msg);
}
