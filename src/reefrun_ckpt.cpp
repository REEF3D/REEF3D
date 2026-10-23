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
#include<csignal>
#include<cmath>
#include<map>
#include<set>
#include<sstream>
#include<algorithm>

using namespace std;

namespace
{
// SIGTERM, SIGUSR1 (a batch system before the end of the job, a container before its restart): the run
// writes a checkpoint at the next step boundary and stops (signal_stop)
volatile sig_atomic_t stop_signal = 0;

extern "C" void on_stop_signal(int s)
{
    stop_signal = s;
}

void stop_signals_ini()
{
    struct sigaction sa;
    sa.sa_handler = on_stop_signal;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = SA_RESTART;
    sigaction(SIGTERM,&sa,nullptr);
    sigaction(SIGUSR1,&sa,nullptr);
}

// the input files of the run, copied into every checkpoint (rank 0): a restart lists what changed
const char *const INPUTS[] = {"ctrl.txt","control.txt"};

void copy_file(const string &from, const string &to)
{
    ifstream i(from.c_str(),ios::binary);
    if(!i)
    return;
    ofstream o(to.c_str(),ios::binary);
    o<<i.rdbuf();
}

// the lines of an input file by key (the letter and the number), comments (/) and blank lines left out
map<string,vector<string>> input_lines(const string &file)
{
    map<string,vector<string>> m;
    ifstream f(file.c_str());
    string l;
    while(getline(f,l))
    {
        const size_t c = l.find('/');
        if(c!=string::npos)
        l.erase(c);
        istringstream is(l);
        string a, b, v, rest;
        if(!(is>>a>>b) || !isalpha((unsigned char)a[0]))
        continue;
        while(is>>v)
        rest += (rest.empty() ? "" : " ")+v;
        m[a+" "+b].push_back(rest);
    }
    for(auto &e : m)
    sort(e.second.begin(),e.second.end());
    return m;
}

// what a change of the setting does on a restart, as far as it was looked at (reefrun_ckpt.h)
string restart_note(const string &key)
{
    const char k = key[0];
    if(k=='P')
    return "output";
    if(key=="M 10")
    return "ranks: any decomposition";
    if(key=="N 41" || key=="N 45")
    return "the end of the run";
    if(key=="N 47" || key=="N 48" || key=="N 49")
    return "the time step: used in practice";
    if(key=="G 45" || key=="G 46" || key=="G 47")
    return "the rebalance";
    if(key=="A 10" || key=="X 10" || key=="Z 40")
    return "the model or its bodies: the restart is refused or does not fit";
    return "not checked for restarts";
}

// rank 0, at a restart: the settings of ctrl.txt and control.txt that differ from those of the checkpoint
void settings_changed(const string &dir)
{
    for(const char *in : INPUTS)
    {
        const string old = dir+"/"+in;
        ifstream t(old.c_str());
        if(!t)
        continue;
        const auto a = input_lines(old), b = input_lines(in);
        vector<string> out;
        set<string> keys;
        for(const auto &e : a) keys.insert(e.first);
        for(const auto &e : b) keys.insert(e.first);
        for(const string &k : keys)
        {
            if(k=="I 45")
            continue;
            const auto ia = a.find(k), ib = b.find(k);
            const vector<string> va = (ia==a.end()) ? vector<string>() : ia->second;
            const vector<string> vb = (ib==b.end()) ? vector<string>() : ib->second;
            if(va==vb)
            continue;
            string s = "  "+k+": ";
            if(va.size()==1 && vb.size()==1)
            s += va[0]+" -> "+vb[0];
            else if(va.empty())
            s += "added"+(vb.size()>1 ? " ("+to_string(vb.size())+" lines)" : (vb[0].empty() ? string() : " "+vb[0]));
            else if(vb.empty())
            s += "removed";
            else
            s += to_string(va.size())+" lines -> "+to_string(vb.size())+" lines";
            out.push_back(s+"   ["+restart_note(k)+"]");
        }
        if(out.empty())
        continue;
        cout<<"I 45: "<<in<<" differs from the one of the checkpoint ("<<out.size()<<(out.size()==1 ? " setting" : " settings")<<"):"<<endl;
        for(const string &s : out)
        cout<<s<<endl;
    }
}

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
    // the pressure reference of the last step (B 30; the loads of the 6DOF bodies take it off)
    f("pressgage",p->pressgage);
    // sediment (S 10): the sediment time and step (S 19, S 29, outputs), the morphological times (S 42),
    // the bed range and the slide cells of the log
    f("sedtime",p->sedtime); f("dtsed",p->dtsed); f("sedsimtime",p->sedsimtime); f("sedwavetime",p->sedwavetime);
    f("bedmin",p->bedmin); f("bedmax",p->bedmax); f("slidecells",p->slidecells);
}
}

reefrun_ckpt::reefrun_ckpt(lexer *pp, ghostcell *ppgc, hooks hh) : p(pp), pgc(ppgc), h(std::move(hh))
{
    // a case the checkpoints do not cover runs on without them: said at the start, not only at the first checkpoint
    const string un = (enabled() && h.unsupported) ? h.unsupported() : string();
    if(p->mpirank==0 && !un.empty())
    cout<<endl<<"!!! P 48: "<<un<<" -- this run writes no checkpoints !!!"<<endl<<endl;

    // SIGTERM, SIGUSR1: a checkpoint at the next step boundary, then the stop (signal_stop); a case not covered
    // keeps the default (the signal ends the run)
    if(enabled() && un.empty())
    stop_signals_ini();
}

bool reefrun_ckpt::enabled() const
{
    return p->P48>0 || p->P48_dt>0.0;
}

int reefrun_ckpt::keep() const
{
    if(p->P49>=0)
    return p->P49;
    return h.keep_last ? 1 : 2;
}

// P 48 N dt: the next multiple of dt after the simulated time
double reefrun_ckpt::next_time() const
{
    const double d = p->P48_dt;
    return (std::floor(p->simtime/d + 1.0e-9) + 1.0)*d;
}

void reefrun_ckpt::step()
{
    bool due = (p->P48>0 && p->count%p->P48==0);

    // P 48 N dt: the first step at or past each multiple of dt (the spacing in simulated time does not
    // change with the time step)
    if(p->P48_dt>0.0)
    {
        if(next_t<=0.0)
        next_t = next_time();
        if(p->simtime >= next_t - 1.0e-9*p->P48_dt)
        {
            due = true;
            next_t = next_time();
        }
    }

    if(due)
    write();
    if(enabled())
    signal_stop();
    if(p->G46>0.0)
    rebalance_check();
}

void reefrun_ckpt::finish()
{
    if(!enabled() || (last_count==p->count && last_sub==sub))
    return;
    write();
}

bool reefrun_ckpt::signal_stop()
{
    int s = stop_signal, all = 0;
    MPI_Allreduce(&s,&all,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    if(all==0)
    return false;

    if(p->mpirank==0)
    cout<<endl<<"signal "<<(all==SIGTERM ? "SIGTERM" : (all==SIGUSR1 ? "SIGUSR1" : to_string(all)))<<": checkpoint at step "<<p->count
        <<", simtime "<<setprecision(10)<<p->simtime<<" and stop (I 45 1 continues)"<<endl;
    if(last_count!=p->count || last_sub!=sub)
    write();
    if(h.flush)
    h.flush();
    MPI_Barrier(MPI_COMM_WORLD);
    pgc->final();
    return true;
}

void reefrun_ckpt::write()
{
    const string un = h.unsupported ? h.unsupported() : string();
    if(!un.empty())
    {
        if(p->mpirank==0 && !told)
        cout<<"!!! P 48: "<<un<<" -- no checkpoint !!!"<<endl;
        told = true;
        return;
    }

    const double t0 = pgc->timer();
    reefckpt ck;
    ck.write_open(p,p->count,sub);

    ck.scalar("count",double(p->count));
    ck.scalar("rebal_L",rebal_L);
    ck.scalar("rebal_n",double(rebal_n));
    ck.scalar("printcount",double(p->printcount));
    ck.scalar("printcount_sixdof",double(p->printcount_sixdof));
    ck.scalar("sediter",double(p->sediter));
    run_scalars(p,[&](const char *n, double &v){ ck.scalar(n,v); });

    vector<double> pv;
    if(h.printer_state)
    h.printer_state(pv);
    ck.blob("printer",pv);

    h.write(ck);

    // the input files: a restart lists what changed since
    if(p->mpirank==0)
    for(const char *in : INPUTS)
    copy_file(in,ck.folder()+"/"+in);

    // the outputs on disk up to the checkpoint: their sizes go with it (a restart cuts them back to these)
    if(h.flush)
    h.flush();

    ck.write_close(p);
    write_s = pgc->timer()-t0;
    last_count = p->count;
    last_sub = sub;

    // P 49: only the newest checkpoints are kept (REEF3D_Checkpoint/history), the older ones are removed once
    // this one is the latest; not given: two (one an instability may not have reached yet), keep_last: one;
    // 0: all (the history is kept all the same)
    if(p->mpirank==0)
    reefckpt::prune(ck.folder(),keep());

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

    // a failed restart of a rebalance (G 46) removes the marker, so that tools/reef3d_rebalance.sh does not start it again
    auto fail = [&](const string &m)
    {
        if(p->mpirank==0)
        {
            cout<<endl<<"!!! I 45: "<<m<<" !!!"<<endl<<endl;
            remove("REEF3D_Checkpoint/rebalance");
        }
        pgc->final();
        exit(1);
    };

    if(!msg.empty() || !ck.read_open(p,"",msg))
    fail(msg);

    if(!h.read(ck,msg))
    fail(msg.empty() ? "checkpoint "+ck.folder()+" incomplete for this case" : msg);

    p->count = int(ck.get_scalar("count"));
    next_t = 0.0;
    last_count = p->count;
    last_sub = sub;
    rebal_L = ck.get_scalar("rebal_L");
    rebal_n = int(ck.get_scalar("rebal_n"));
    p->count_statestart = p->count;
    p->printcount = int(ck.get_scalar("printcount"));
    p->printcount_sixdof = int(ck.get_scalar("printcount_sixdof"));
    p->sediter = int(ck.get_scalar("sediter",p->sediter));
    run_scalars(p,[&](const char *n, double &v){ v = ck.get_scalar(n,v); });
    if(h.printer_restore)
    h.printer_restore(ck.get_blob("printer"));

    // a rebalance restart (G 46): the marker is used up
    if(p->mpirank==0)
    remove("REEF3D_Checkpoint/rebalance");

    // the checkpoint of the restart is the newest of the history (P 49): the ones after it (a restart from an
    // earlier checkpoint, I 45 <count>) are removed first
    if(p->mpirank==0)
    reefckpt::prune(ck.folder(),0);

    if(p->mpirank==0)
    settings_changed(ck.folder());

    if(p->mpirank==0)
    cout<<"restart from checkpoint "<<ck.folder()<<" (I 45): count "<<p->count<<", simtime "<<setprecision(10)<<p->simtime
        <<(ck.same_decomposition() ? ", same decomposition" : ", decomposition changed from "+to_string(ck.pieces())+" pieces")<<endl;
}

// G 46: every G 47 steps the load of the ranks since the last check (busy max/mean, the MPI waits
// taken out).  Above G 46 the run stops to rebalance: the load map of the mesh refinement, a
// checkpoint, the marker REEF3D_Checkpoint/rebalance, exit code 3.  The next start (the same command,
// tools/reef3d_rebalance.sh repeats it) sees the marker: I 45 1 with the weighted partition of the
// load map (M 41 1).  If the first check after a rebalance restart is not 10 % better than the load
// that stopped the run, there is no further rebalance in this run.  A module that knows the load of a
// partition from its load map (hooks.predict: SEASTATE) stops only if the weighted partition is
// 10 % better than the present one by that map, and (hooks.pays) only if the time it saves over the
// rest of the run exceeds the cost of the stop, the checkpoint and the new start.
void reefrun_ckpt::rebalance_check()
{
    const string un = h.unsupported ? h.unsupported() : string();
    if(!un.empty())
    {
        if(p->mpirank==0 && p->count==p->G47)
        cout<<"G 46: "<<un<<" -- no rebalance"<<endl;
        return;
    }

    rebalance(p->G47>0 && p->count%p->G47==0, p->count+p->G47>p->N45 || (p->dt>0.0 && p->simtime+p->G47*p->dt>p->N41));
}

void reefrun_ckpt::rebalance_check_iter(int k, int kmax)
{
    const string un = h.unsupported ? h.unsupported() : string();
    if(!un.empty())
    {
        if(p->mpirank==0 && k==p->G47)
        cout<<"G 46: "<<un<<" -- no rebalance"<<endl;
        return;
    }

    rebalance(p->G47>0 && k%p->G47==0, k+p->G47>kmax);
}

// due: the check is due now; last: the run (or the iterations of the step) end within G 47
void reefrun_ckpt::rebalance(bool due, bool last)
{
    if(p->mpi_size<2 || rebal_off || !due)
    return;

    const double L = h.imbalance ? h.imbalance() : amr_window_imbalance();
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
    if(last)
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
    // 1: it differs; 2: it differs, but would not lower the load enough (predict); 3: the restart would not pay (pays)
    int differs = 0;
    double Lnew = 0.0;
    if(p->mpirank==0)
    {
        vector<vector<char>> rb;
        vector<int> nb;
        string msg;
        if(dmpart::build("DIVEMesh_Grid/grid-global.dat",p->mpi_size,rb,msg,1,&nb,true))
        {
            differs = (nb!=now) ? 1 : 0;
            if(differs && h.predict)
            {
                Lnew = h.predict(nb);
                if(Lnew>0.9*L)
                differs = 2;
            }
            // 3: the restart would not pay
            if(differs==1 && h.pays)
            {
                string why;
                if(!h.pays(L,Lnew,why))
                differs = 3;
                if(!why.empty())
                cout<<"G 46: "<<why<<endl;
            }
        }
        else
        cout<<"G 46: no rebalance -- "<<msg<<" (needs DIVEMesh M 50 1, 2)"<<endl;
    }
    MPI_Bcast(&differs,1,MPI_INT,0,MPI_COMM_WORLD);
    if(differs!=1)
    {
        if(p->mpirank==0 && differs==0)
        cout<<"G 46: the weighted partition of the load map is the present one -- no rebalance"<<endl;
        if(p->mpirank==0 && differs==3)
        cout<<"G 46: the restart does not pay -- no rebalance"<<endl;
        if(p->mpirank==0 && differs==2)
        {
            const streamsize pr = cout.precision();
            cout<<"G 46: the weighted partition of the load map would give the load max/mean "<<setprecision(3)<<Lnew
                <<", less than 10 % below the present one -- no rebalance"<<endl;
            cout.precision(pr);
        }
        return;
    }
    if(p->mpirank==0 && h.predict)
    {
        const streamsize pr = cout.precision();
        cout<<"G 46: the weighted partition of the load map gives the load max/mean "<<setprecision(3)<<Lnew<<" by the load map"<<endl;
        cout.precision(pr);
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
