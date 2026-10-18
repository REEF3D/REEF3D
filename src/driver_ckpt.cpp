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

//  Checkpoints on any decomposition (P 48, I 45) and the rebalance (G 46) of FNPF and NHFLOW: the
//  state of the modules (the run part: reefrun_ckpt.h).  The arrays of the model (fnpf_ckpt_fields.h,
//  nhflow_ckpt_fields.h, the momentum and turbulence objects), the patches of the mesh refinement and
//  the counters of the outputs.  On the same decomposition the restart continues bitwise; on another
//  one the ghost cells come from the neighbours' interiors and the domain boundaries of the writing run.
//  Not there yet (refused with a message): 6DOF bodies; NHFLOW with sediment, Lagrangian particles,
//  DEM, FEM or CORAL structures, SEASTATE, membranes, multi-block grids.

#include"driver.h"
#include"lexer.h"
#include"ghostcell.h"
#include"fdm_fnpf.h"
#include"printer.h"
#include"fnpf_amr.h"
#include"reefckpt.h"
#include"reefrun_ckpt.h"
#include"fnpf_ckpt_fields.h"
#include"nhflow_ckpt_fields.h"
#include"fdm_nhf.h"
#include"nhflow_amr.h"
#include"nhflow_turbulence.h"
#include"nhflow_momentum.h"
#include"ioflow.h"

using namespace std;

// what the checkpoints do not cover yet, empty: all
string driver::ckpt_unsupported()
{
    if(p->A10!=3 && p->A10!=5)
    return "checkpoints are there for SFLOW, FNPF and NHFLOW only";
    if(p->X10>0)
    return "6DOF bodies are not in the checkpoints yet";
    if(p->A10==5)
    {
        if(p->S10>0) return "NHFLOW sediment is not in the checkpoints yet";
        if(p->L10>0) return "NHFLOW Lagrangian particles are not in the checkpoints yet";
        if(p->E10>0) return "DEM is not in the checkpoints yet";
        if(p->Z30>0 || p->Z40>0) return "FEM and CORAL structures are not in the checkpoints yet";
        if(p->A750>0) return "SEASTATE is not in the checkpoints yet";
        if(p->X330>0) return "membranes are not in the checkpoints yet";
        if(p->geo_mb>0) return "multi-block grids are not in the checkpoints yet";
    }
    return "";
}

reefrun_ckpt& driver::ckpt()
{
    if(rck!=nullptr)
    return *rck;

    reefrun_ckpt::hooks h;

    h.unsupported = [this](){ return ckpt_unsupported(); };

    h.write = [this](reefckpt &ck)
    {
        auto wd = [&](const char *n, char k, double *V){ if(V==nullptr) return; if(k=='S') ck.put_S(n,V); else if(k=='F') ck.put_F(n,V); else ck.put_C(n,V); };
        auto wi = [&](const char *n, char k, int *V)   { if(V==nullptr) return; if(k=='S') ck.put_Si(n,V); else if(k=='F') ck.put_Fi(n,V); else ck.put_Ci(n,V); };

        if(p->A10==3)
        {
            fnpf_ckpt_fields(p,c,wd,wi);
            if(pfamr!=nullptr)
            pfamr->ckpt_write(ck);
        }
        if(p->A10==5)
        {
            nhflow_ckpt_fields(p,d,wd,wi);
            pnhfturb->ckpt_fields(wd,wi);
            pnhfmom->ckpt_fields(wd,wi);
            if(pnhfamr!=nullptr)
            pnhfamr->ckpt_write(ck);
        }

        // the state of the boundary conditions (local to the rank)
        vector<double> io;
        pflow->ckpt_state(p,io);
        ck.rank_blob("ioflow",io);
    };

    h.read = [this](reefckpt &ck, string &msg)
    {
        bool ok = true;
        auto rd = [&](const char *n, char k, double *V){ if(V==nullptr || !ck.has(n)) return; ok = (k=='S' ? ck.get_S(n,V) : (k=='F' ? ck.get_F(n,V) : ck.get_C(n,V))) && ok; };
        auto ri = [&](const char *n, char k, int *V)   { if(V==nullptr || !ck.has(n)) return; ok = (k=='S' ? ck.get_Si(n,V) : (k=='F' ? ck.get_Fi(n,V) : ck.get_Ci(n,V))) && ok; };
        auto level0 = [&]()
        {
            if(p->A10==3)
            fnpf_ckpt_fields(p,c,rd,ri);
            if(p->A10==5)
            {
                nhflow_ckpt_fields(p,d,rd,ri);
                pnhfturb->ckpt_fields(rd,ri);
                pnhfmom->ckpt_fields(rd,ri);
            }
        };

        level0();
        if(p->A10==5)
        pnhfmom->ckpt_restored();
        if(p->A10==3)
        ok = ck.has("eta") && ck.has("Fifsf") && ck.has("Fi") && ok;
        if(p->A10==5)
        ok = ck.has("eta") && ck.has("UH") && ck.has("P") && ok;

        // mesh refinement: the patches, then level 0 once more (the regrid to the layout of the
        // checkpoint restricted its patches)
        if(pfamr!=nullptr || pnhfamr!=nullptr)
        {
            const bool aok = (pfamr!=nullptr) ? pfamr->ckpt_read(ck,msg) : pnhfamr->ckpt_read(ck,msg);
            if(!aok)
            return false;
            level0();
        }

        // the state of the boundary conditions: on the same decomposition (else they start again)
        if(ck.same_decomposition())
        pflow->ckpt_restore(p,ck.get_rank_blob("ioflow",p->mpirank));
        return ok;
    };

    h.printer_state = [this](vector<double> &v){ pprint->ckpt_state(v); };
    h.printer_restore = [this](const vector<double> &v){ pprint->ckpt_restore(v); };

    h.flush = [this]()
    {
        pprint->flush_outputs();
        if(pfamr!=nullptr)
        pfamr->flush_outputs();
        if(pnhfamr!=nullptr)
        pnhfamr->flush_outputs();
    };

    h.loadmap = [this](lexer *pp)
    {
        if(pfamr!=nullptr)
        pfamr->write_loadmap(pp);
        if(pnhfamr!=nullptr)
        pnhfamr->write_loadmap(pp);
    };

    rck = new reefrun_ckpt(p,pgc,std::move(h));
    return *rck;
}

void driver::ckpt_write()
{
    ckpt().write();
}

void driver::ckpt_read()
{
    ckpt().read();
}

void driver::rebalance_check()
{
    ckpt().rebalance_check();
}
