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

//  FNPF checkpoints on any decomposition (P 48: write every P 48 steps, I 45 1: restart from the
//  latest; reefckpt.h).  The state at the end of a step: the counters and times of the run, every
//  slice of fdm_fnpf, the potential and velocities, the sigma metrics and the flags of the lexer
//  that change in time (wetting and drying), and the counters of the outputs.  On the same
//  decomposition the restart continues bitwise; on another one the ghost cells come from the
//  neighbours' interiors and the domain boundaries of the writing run.

#include"driver.h"
#include"lexer.h"
#include"ghostcell.h"
#include"fdm_fnpf.h"
#include"printer.h"
#include"fnpf_amr.h"
#include"reefckpt.h"
#include"fnpf_ckpt_fields.h"
#include<iostream>
#include<iomanip>
#include<fstream>
#include<cstdio>
#include<mpi.h>
#include"amr_load.h"
#include"dmpart_api.h"

using namespace std;

namespace
{
// the scalars of the run
template<class F>
void run_scalars(lexer *p, F &&f)
{
    f("simtime",p->simtime); f("dt",p->dt); f("dt_old",p->dt_old);
    f("printtime",p->printtime); f("sedprinttime",p->sedprinttime); f("fsfprinttime",p->fsfprinttime);
    f("fsfsedprinttime",p->fsfsedprinttime); f("probeprinttime",p->probeprinttime);
    f("stateprinttime",p->stateprinttime); f("exportprinttime",p->exportprinttime); f("wavetime",p->wavetime);
    f("totaltime",p->totaltime); f("gctotaltime",p->gctotaltime); f("Xtotaltime",p->Xtotaltime);
}
}

void driver::fnpf_ckpt_write()
{
    if(p->X10>0)
    {
        if(p->mpirank==0 && p->count==p->P48)
        cout<<"P 48: FNPF checkpoints with 6DOF bodies are not there yet -- no checkpoint"<<endl;
        return;
    }

    const double t0 = pgc->timer();
    reefckpt ck;
    ck.write_open(p,p->count);

    ck.scalar("count",double(p->count));
    ck.scalar("rebal_L",rebal_L);
    ck.scalar("rebal_n",double(rebal_n));
    ck.scalar("printcount",double(p->printcount));
    ck.scalar("printcount_sixdof",double(p->printcount_sixdof));
    run_scalars(p,[&](const char *n, double &v){ ck.scalar(n,v); });

    vector<double> pv;
    pprint->ckpt_state(pv);
    ck.blob("printer",pv);

    fnpf_ckpt_fields(p,c,
        [&](const char *n, char k, double *V){ if(V==nullptr) return; if(k=='S') ck.put_S(n,V); else if(k=='F') ck.put_F(n,V); else ck.put_C(n,V); },
        [&](const char *n, char k, int *V)   { if(V==nullptr) return; if(k=='S') ck.put_Si(n,V); else if(k=='F') ck.put_Fi(n,V); else ck.put_Ci(n,V); });

    // mesh refinement: the patches (fnpf_amr_ckpt.cpp)
    if(pfamr!=nullptr)
    pfamr->ckpt_write(ck);

    ck.write_close(p);

    if(p->mpirank==0)
    {
        const streamsize pr = cout.precision();
        cout<<"checkpoint "<<ck.folder()<<" (P 48) "<<setprecision(3)<<pgc->timer()-t0<<" s"<<endl;
        cout.precision(pr);
    }
}

void driver::fnpf_ckpt_read()
{
    reefckpt ck;
    string msg;

    if(p->X10>0)
    msg = "FNPF restart with 6DOF bodies is not there yet";

    if(!msg.empty() || !ck.read_open(p,"",msg))
    {
        if(p->mpirank==0)
        cout<<endl<<"!!! I 45 1: "<<msg<<" !!!"<<endl<<endl;
        pgc->final();
        exit(1);
    }

    bool ok = true;
    fnpf_ckpt_fields(p,c,
        [&](const char *n, char k, double *V){ if(V==nullptr || !ck.has(n)) return; ok = (k=='S' ? ck.get_S(n,V) : (k=='F' ? ck.get_F(n,V) : ck.get_C(n,V))) && ok; },
        [&](const char *n, char k, int *V)   { if(V==nullptr || !ck.has(n)) return; ok = (k=='S' ? ck.get_Si(n,V) : (k=='F' ? ck.get_Fi(n,V) : ck.get_Ci(n,V))) && ok; });
    ok = ck.has("eta") && ck.has("Fifsf") && ck.has("Fi") && ok;

    // mesh refinement: the patches
    if(pfamr!=nullptr)
    {
        string amsg;
        if(!pfamr->ckpt_read(ck,amsg))
        {
            if(p->mpirank==0)
            cout<<endl<<"!!! I 45 1: "<<amsg<<" !!!"<<endl<<endl;
            pgc->final();
            exit(1);
        }
    }

    if(!ok)
    {
        if(p->mpirank==0)
        cout<<endl<<"!!! I 45 1: checkpoint "<<ck.folder()<<" incomplete for this case !!!"<<endl<<endl;
        pgc->final();
        exit(1);
    }

    p->count = int(ck.get_scalar("count"));
    rebal_L = ck.get_scalar("rebal_L");
    rebal_n = int(ck.get_scalar("rebal_n"));
    p->count_statestart = p->count;
    p->printcount = int(ck.get_scalar("printcount"));
    p->printcount_sixdof = int(ck.get_scalar("printcount_sixdof"));
    run_scalars(p,[&](const char *n, double &v){ v = ck.get_scalar(n,v); });
    pprint->ckpt_restore(ck.get_blob("printer"));

    // a rebalance restart (G 46): the marker is used up
    if(p->mpirank==0)
    remove("REEF3D_Checkpoint/rebalance");

    if(p->mpirank==0)
    cout<<"restart from checkpoint "<<ck.folder()<<" (I 45): count "<<p->count<<", simtime "<<setprecision(10)<<p->simtime
        <<(ck.same_decomposition() ? ", same decomposition" : ", decomposition changed from "+to_string(ck.pieces())+" pieces")<<endl;
}

// G 46: every G 47 steps the load of the ranks since the last check (busy max/mean, the MPI waits
// taken out).  Above G 46 the run stops to rebalance: the load map of the mesh refinement, a
// checkpoint, the marker REEF3D_Checkpoint/rebalance, exit code 3.  The next start (the same command,
// tools/reef3d_rebalance.sh repeats it) sees the marker: I 45 1 with the weighted partition of the
// load map (M 41 1).  If the first check after a rebalance restart is not 10 % better than the load
// that stopped the run, there is no further rebalance in this run.
void driver::fnpf_rebalance_check()
{
    if(p->mpi_size<2 || rebal_off || p->G47<=0 || p->count%p->G47!=0)
    return;

    const double L = amr_window_imbalance();
    if(L<=0.0)
    return;

    if(p->mpirank==0)
    {
        const streamsize pr = cout.precision();
        cout<<"G 46: load of the ranks max/mean "<<setprecision(3)<<L<<" (threshold "<<p->G46<<")"<<endl;
        cout.precision(pr);
    }

    // the first check after a rebalance restart
    if(rebal_L>0.0)
    {
        if(L>0.9*rebal_L)
        {
            rebal_off = true;
            if(p->mpirank==0)
            cout<<"G 46: the rebalance did not lower the load max/mean enough ("<<rebal_L<<" -> "<<L<<") -- no further rebalance in this run"<<endl;
        }
        rebal_L = 0.0;
        return;
    }

    if(L<=p->G46)
    return;

    // the load map, and the weighted partition it gives (rank 0, in place): a stop only if that
    // decomposition differs from the present one
    if(pfamr!=nullptr)
    pfamr->write_loadmap(p);

    int mybox[6] = {p->origin_i,p->origin_j,p->origin_k,p->knox,p->knoy,p->knoz};
    vector<int> now(6*size_t(p->mpi_size));
    MPI_Gather(mybox,6,MPI_INT,now.data(),6,MPI_INT,0,MPI_COMM_WORLD);
    int differs = 0;
    if(p->mpirank==0)
    {
        vector<vector<char>> rb;
        vector<int> nb;
        string msg;
        if(dmpart::build("DIVEMesh_Grid/grid-global.dat",p->mpi_size,rb,msg,1,&nb,true))
        differs = (nb!=now) ? 1 : 0;
        else
        cout<<"G 46: no rebalance -- "<<msg<<" (needs DIVEMesh M 50 1, 2)"<<endl;
    }
    MPI_Bcast(&differs,1,MPI_INT,0,MPI_COMM_WORLD);
    if(!differs)
    {
        if(p->mpirank==0)
        cout<<"G 46: the weighted partition of the load map is the present one -- no rebalance"<<endl;
        return;
    }

    if(p->mpirank==0)
    cout<<"G 46: rebalance -- checkpoint and stop (exit code 3); the next start decomposes for the measured load"<<endl;

    rebal_L = L;
    ++rebal_n;
    fnpf_ckpt_write();

    pprint->flush_outputs();
    if(pfamr!=nullptr)
    pfamr->flush_outputs();

    if(p->mpirank==0)
    {
        ofstream m("REEF3D_Checkpoint/rebalance");
        m<<p->count<<" "<<L<<endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    pgc->final();
    exit(3);
}
