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

//  SFLOW checkpoints on any decomposition (P 48, I 45) and the rebalance (G 46): the state of the
//  module for the run part (reefrun_ckpt.h).  The slices of fdm2D and the wet/dry flags
//  (sflow_ckpt_fields.h), the arrays of the momentum, turbulence, pressure, diffusion and free-surface
//  objects, the patches of the mesh refinement and the counters of the outputs.  On the same
//  decomposition the restart continues bitwise.
//  The blocks 1.. of a multi-block grid are patches of the mesh refinement (reefamr_ckpt.cpp).
//  The 6DOF bodies (one-way bodies X 10 2, ship waves X 10 3): their state as a blob, the pressure patch
//  of the ship waves as a slice.
//  Not there yet (refused with a message): 6DOF wavemaker objects, sediment, SEASTATE.

#include"sflow_f.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"reefckpt.h"
#include"reefrun_ckpt.h"
#include"sflow_ckpt_fields.h"
#include"sflow_momentum.h"
#include"sflow_turbulence.h"
#include"sflow_pressure.h"
#include"sflow_diffusion.h"
#include"sflow_fsf.h"
#include"sflow_vtp_fsf.h"
#include"sflow_vtp_bed.h"
#include"sflow_amr.h"
#include"ioflow.h"
#include"6DOF.h"
#include"cfd_ckpt.h"
#include"sediment.h"

using namespace std;

string sflow_f::ckpt_unsupported(lexer *p)
{
    if(p->X10==4)
    return "6DOF wavemaker objects (X 10 4) are not in the checkpoints yet";
    if(p->S10>0 && (p->G1>0 || p->geo_mb>0))
    return "SFLOW sediment with mesh refinement or multi-block grids is not in the checkpoints yet";
    if(p->A750>0)
    return "SEASTATE is not in the checkpoints yet";
    return "";
}

reefrun_ckpt& sflow_f::ckpt(lexer *p, fdm2D *b, ghostcell *pgc)
{
    if(rck!=nullptr)
    return *rck;

    reefrun_ckpt::hooks h;

    h.unsupported = [this,p](){ return ckpt_unsupported(p); };

    auto fields = [this,p,b](auto &fd, auto &fi)
    {
        sflow_ckpt_fields(p,b,fd,fi);
        pmom->ckpt_fields(fd,fi);
        pturb->ckpt_fields(fd,fi);
        ppress->ckpt_fields(fd,fi);
        pdiff->ckpt_fields(fd,fi);
        pfsf->ckpt_fields(fd,fi);
        if(pamr!=nullptr)
        pamr->ckpt_fields0(fd,fi);
        // the 6DOF bodies on the grid
        if(p->X10>0)
        p6dof->ckpt_fields(fd,fi);
    };

    h.write = [this,p,fields](reefckpt &ck)
    {
        std::function<void(const char*,char,double*)> wd = [&](const char *n, char, double *V){ if(V!=nullptr) ck.put_S(n,V); };
        std::function<void(const char*,char,int*)> wi = [&](const char *n, char, int *V){ if(V!=nullptr) ck.put_Si(n,V); };
        fields(wd,wi);
        if(pamr!=nullptr)
        pamr->ckpt_write(ck);

        // the 6DOF bodies (the same on all ranks) and the pressure patch of the ship waves
        if(p->X10>0)
        {
            vector<double> sd;
            p6dof->ckpt_state(p,sd);
            ck.blob("sixdof",sd);
            if(p6dof->ckpt_slice()!=nullptr)
            ck.put_S("sixdof_slice",p6dof->ckpt_slice());
        }

        // sediment (S 10): the registry of the grid (the slices of the sediment objects), the scalars, those
        // of the rank, the sediment cell flags
        if(p->S10>0)
        {
            cfd_ckpt_write(p,ck,"reg");
            vector<double> sv, sr;
            psed->ckpt_state(p,sv);
            ck.blob("sediment",sv);
            psed->ckpt_rank_state(p,sr);
            ck.rank_blob("sediment_rank",sr);
            if(p->DFBED!=nullptr)
            ck.put_Si("DFBED",p->DFBED);
        }

        // the state of the boundary conditions (local to the rank)
        vector<double> io;
        pflow->ckpt_state(p,io);
        ck.rank_blob("ioflow",io);
    };

    h.read = [this,p,pgc,fields](reefckpt &ck, string &msg)
    {
        bool ok = true;
        std::function<void(const char*,char,double*)> rd = [&](const char *n, char, double *V){ if(V==nullptr || !ck.has(n)) return; ok = ck.get_S(n,V) && ok; };
        std::function<void(const char*,char,int*)> ri = [&](const char *n, char, int *V){ if(V==nullptr || !ck.has(n)) return; ok = ck.get_Si(n,V) && ok; };

        // sediment: the registry of the grid first, the explicit lists after it
        if(p->S10>0 && !cfd_ckpt_read(p,ck,msg,"reg"))
        return false;

        fields(rd,ri);
        ok = ck.has("eta") && ck.has("UH") && ck.has("WL") && ok;
        pmom->ckpt_restored();

        // mesh refinement: the patches, then level 0 once more (the regrid to the layout of the
        // checkpoint restricted its patches)
        if(pamr!=nullptr)
        {
            if(!pamr->ckpt_read(ck,msg))
            return false;
            fields(rd,ri);
        }

        // the 6DOF bodies: their state, the trimesh at the pose
        if(p->X10>0)
        {
            if(!ck.has("sixdof") || !p6dof->ckpt_restore(p,pgc,ck.get_blob("sixdof")))
            {
                msg = "the 6DOF bodies of the checkpoint do not fit this run";
                return false;
            }
            if(p6dof->ckpt_slice()!=nullptr && ck.has("sixdof_slice"))
            ok = ck.get_S("sixdof_slice",p6dof->ckpt_slice()) && ok;
        }

        // sediment (S 10)
        if(p->S10>0)
        {
            if(!ck.has("sediment") || !psed->ckpt_restore(p,pgc,ck.get_blob("sediment")))
            {
                msg = "the checkpoint has no sediment state that fits this run (S 10)";
                return false;
            }
            const vector<vector<double>> sr = ck.same_decomposition() ? vector<vector<double>>{ck.get_rank_blob("sediment_rank",p->mpirank)}
                                                                      : ck.get_rank_blobs("sediment_rank");
            psed->ckpt_rank_restore(p,pgc,sr,ck.same_decomposition());
            if(p->DFBED!=nullptr && ck.has("DFBED"))
            ok = ck.get_Si("DFBED",p->DFBED) && ok;
        }

        // the state of the boundary conditions: on the same decomposition (else they start again)
        if(ck.same_decomposition())
        pflow->ckpt_restore(p,ck.get_rank_blob("ioflow",p->mpirank));
        return ok;
    };

    h.printer_state = [this](vector<double> &v){ pprintbed->ckpt_state(v); };
    h.printer_restore = [this](const vector<double> &v){ size_t k=0; pprintbed->ckpt_restore(v,k); };

    h.flush = [this]()
    {
        pprint->flush_outputs();
        if(pamr!=nullptr)
        pamr->flush_outputs();
        mainlogout.flush();
    };

    h.loadmap = [this](lexer *pp)
    {
        if(pamr!=nullptr)
        pamr->write_loadmap(pp);
    };

    rck = new reefrun_ckpt(p,pgc,std::move(h));
    return *rck;
}
