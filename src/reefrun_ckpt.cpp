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

#include"reefrun_ckpt.h"
#include"reefckpt.h"
#include"lexer.h"
#include"ghostcell.h"
#include"amr_load.h"
#include"dmpart_api.h"
#include<iostream>
#include<iomanip>
#include<fstream>
#include<cstdio>
#include<mpi.h>

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
    // the in- and outflow of the last step (ioflow: the turbulence inflow profile uses Ui before the
    // momentum step sets it again), the maxima of the last step
    f("Ui",p->Ui); f("Ua",p->Ua); f("Uo",p->Uo); f("Hi",p->Hi); f("Ho",p->Ho); f("Qi",p->Qi); f("Qo",p->Qo);
    f("umax",p->umax); f("vmax",p->vmax); f("wmax",p->wmax); f("viscmax",p->viscmax);
    f("kinmax",p->kinmax); f("epsmax",p->epsmax); f("omegamax",p->omegamax);
    f("pressmin",p->pressmin); f("pressmax",p->pressmax); f("mindt",p->mindt); f("maxdt",p->maxdt);
    f("phimean",p->phimean); f("phiout",p->phiout); f("phiin",p->phiin);
}
}

reefrun_ckpt::reefrun_ckpt(lexer *pp, ghostcell *ppgc, hooks hh) : p(pp), pgc(ppgc), h(std::move(hh))
{
}

void reefrun_ckpt::step()
{
    if(p->P48>0 && p->count%p->P48==0)
    write();
    if(p->G46>0.0)
    rebalance_check();
}

void reefrun_ckpt::write()
{
    const string un = h.unsupported ? h.unsupported() : string();
    if(!un.empty())
    {
        if(p->mpirank==0 && p->count==p->P48)
        cout<<"P 48: "<<un<<" -- no checkpoint"<<endl;
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
    if(h.printer_state)
    h.printer_state(pv);
    ck.blob("printer",pv);

    h.write(ck);

    ck.write_close(p);

    if(p->mpirank==0)
    {
        const streamsize pr = cout.precision();
        cout<<"checkpoint "<<ck.folder()<<" (P 48) "<<setprecision(3)<<pgc->timer()-t0<<" s"<<endl;
        cout.precision(pr);
    }
}

void reefrun_ckpt::read()
{
    reefckpt ck;
    string msg = h.unsupported ? h.unsupported() : string();

    if(!msg.empty() || !ck.read_open(p,"",msg))
    {
        if(p->mpirank==0)
        cout<<endl<<"!!! I 45 1: "<<msg<<" !!!"<<endl<<endl;
        pgc->final();
        exit(1);
    }

    if(!h.read(ck,msg))
    {
        if(p->mpirank==0)
        cout<<endl<<"!!! I 45 1: "<<(msg.empty() ? "checkpoint "+ck.folder()+" incomplete for this case" : msg)<<" !!!"<<endl<<endl;
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
    if(h.printer_restore)
    h.printer_restore(ck.get_blob("printer"));

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
void reefrun_ckpt::rebalance_check()
{
    const string un = h.unsupported ? h.unsupported() : string();
    if(!un.empty())
    {
        if(p->mpirank==0 && p->count==p->G47)
        cout<<"G 46: "<<un<<" -- no rebalance"<<endl;
        return;
    }

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

    // not in the last window of the run (N 45 steps, N 41 time): the restart would not pay
    if(p->count+p->G47>p->N45 || (p->dt>0.0 && p->simtime+p->G47*p->dt>p->N41))
    {
        if(p->mpirank==0)
        cout<<"G 46: the run ends within G 47 steps -- no rebalance"<<endl;
        return;
    }

    // the load map, and the weighted partition it gives (rank 0, in place): a stop only if that
    // decomposition differs from the present one
    if(h.loadmap)
    h.loadmap(p);

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
    write();

    if(h.flush)
    h.flush();

    if(p->mpirank==0)
    {
        ofstream m("REEF3D_Checkpoint/rebalance");
        m<<p->count<<" "<<L<<endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    pgc->final();
    exit(3);
}
