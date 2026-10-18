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
//  Not there yet (refused with a message): 6DOF bodies, sediment, SEASTATE, multi-block grids.

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

using namespace std;

string sflow_f::ckpt_unsupported(lexer *p)
{
    if(p->X10>0)
    return "6DOF bodies are not in the checkpoints yet";
    if(p->S10>0)
    return "SFLOW sediment is not in the checkpoints yet";
    if(p->A750>0)
    return "SEASTATE is not in the checkpoints yet";
    if(p->geo_mb>0)
    return "multi-block grids are not in the checkpoints yet";
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
    };

    h.write = [this,p,fields](reefckpt &ck)
    {
        std::function<void(const char*,char,double*)> wd = [&](const char *n, char, double *V){ if(V!=nullptr) ck.put_S(n,V); };
        std::function<void(const char*,char,int*)> wi = [&](const char *n, char, int *V){ if(V!=nullptr) ck.put_Si(n,V); };
        fields(wd,wi);
        if(pamr!=nullptr)
        pamr->ckpt_write(ck);

        // the state of the boundary conditions (local to the rank)
        vector<double> io;
        pflow->ckpt_state(p,io);
        ck.rank_blob("ioflow",io);
    };

    h.read = [this,p,fields](reefckpt &ck, string &msg)
    {
        bool ok = true;
        std::function<void(const char*,char,double*)> rd = [&](const char *n, char, double *V){ if(V==nullptr || !ck.has(n)) return; ok = ck.get_S(n,V) && ok; };
        std::function<void(const char*,char,int*)> ri = [&](const char *n, char, int *V){ if(V==nullptr || !ck.has(n)) return; ok = ck.get_Si(n,V) && ok; };

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
