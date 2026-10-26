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
Architect: Hans Bihs
--------------------------------------------------------------------*/

#include"seastate_f.h"
#include"fdm_seastate.h"
#include"seastate_grid.h"
#include"seastate_store.h"
#include"seastate_vtp.h"
#include"seastate_exchange.h"
#include"seastate_implicit.h"
#include"seastate_source.h"
#include"seastate_surfbeat.h"
#include"seastate_roller.h"
#include"seastate_forcing.h"
#include"seastate_amr.h"
#include"seastate_bathy.h"
#include"slice4.h"
#include"regression_dump.h"
#include"lexer.h"
#include"ghostcell.h"
#include"runlog.h"
#include"seastate_fas2.h"
#include"reefrun_ckpt.h"
#include<iostream>
#include<mpi.h>
#include<sys/stat.h>
#include<iomanip>

seastate_f::seastate_f(lexer *p, ghostcell *pgc) : pprint(nullptr), preg(nullptr), pex(nullptr), psolv(nullptr), N0(nullptr), psrc(nullptr),
                                                  iter_max(1), iter_done(0), dtw(0.0), coupled(false), conv(0.0),
                                                  etot(0.0), hsmax(0.0), hsmean(0.0), nmin(0.0), cells_active(0.0),
                                                  starttime(0.0), endtime(0.0),
                                                  sb(nullptr), gin(nullptr), proll(nullptr), tnew(0.0), trep(0.0),
                                                  bser(nullptr), wmet(nullptr), tref(0.0), wU10(nullptr), wdir(nullptr), wUx(nullptr), wUy(nullptr),
                                                  pamr(nullptr), bathy(nullptr)
{
    side[0]=side[1]=side[2]=side[3]=0;

    e = new fdm_seastate(p);
}

seastate_f::~seastate_f()
{
    delete rck;
    delete wk0;
    delete pfas;
    delete pamr;
    delete bathy;
    delete pprint;
    delete preg;
    delete pex;
    delete psolv;
    delete N0;
    delete psrc;
    delete proll;
    delete bser;
    delete wU10;
    delete wdir;
    delete wUx;
    delete wUy;
    delete sb;
    delete gin;
    delete e->cg;
    delete e->kw;
    delete e->N;
    delete e->grid;
    delete e;
}

void seastate_f::start(lexer *p, ghostcell *pgc)
{
    const double t0 = pgc->timer();
    ini(p,pgc);

    // Phase 18: the checkpoints and the rebalance (P 48, G 46, G 45); I 45 1: from the latest checkpoint, a stationary
    // step continues at its iteration (the main loop counts it again)
    if(p->P48>0 || p->P48_dt>0.0 || p->G46>0.0 || restart)
    ckpt(p,pgc);
    if(restart)
    {
    rck->read();
    if(ss_resume>0)
    --p->count;
    }
    t_ini = pgc->timer()-t0;

    preg = new regression_dump(p);
    preg->seastate_ini(p,e,pgc);

    if(p->mpirank==0)
    cout<<"starting mainloop.SEASTATE"<<endl;

//-----------MAINLOOP SEASTATE----------------------------
    while(p->count<p->N45 && p->simtime<p->N41)
    {
        ++p->count;
        starttime=pgc->timer();

        if(p->mpirank==0 && (p->count%p->P12==0))
        {
        cout<<"------------------------------------"<<endl;
        cout<<p->count<<endl;
        cout<<"simtime: "<<p->simtime<<endl;
        cout<<"timestep: "<<p->dt<<endl;
        }

        step(p,pgc);

        p->simtime+=p->dt;

        // 2D spectra for the wave generation of FNPF/NHFLOW
        if(p->A760>0)
        handover(p,pgc);

        // printer (mesh refinement: the patches with every print of level 0)
        double ptime=pgc->timer();
        const int pc = p->printcount;
        pprint->start(p,e,pgc);
        if(pamr!=nullptr && p->printcount!=pc)
        pamr->print(p,pgc);
        if(pamr!=nullptr)
        pamr->log(p,pgc,iter_done);
        log_step(p);
        preg->seastate_step(p,e,pgc);
        p->printouttime=pgc->timer()-ptime;

        endtime=pgc->timer();
        p->itertime=endtime-starttime;
        p->totaltime+=p->itertime;
        p->meantime=p->totaltime/double(p->count);

        if(p->mpirank==0 && (p->count%p->P12==0))
        {
        cout<<"Hs max: "<<setprecision(5)<<hsmax<<"   Hs mean: "<<setprecision(5)<<hsmean<<"   iterations: "<<iter_done;
        if(p->A700==2)
        cout<<"   max. relative change of Hs: "<<setprecision(3)<<conv;
        cout<<endl;
        cout<<"printouttime: "<<setprecision(3)<<p->printouttime<<endl;
        cout<<"total time: "<<setprecision(6)<<p->totaltime<<"   average time: "<<setprecision(3)<<p->meantime<<endl;
        }

        p->gctime=0.0;
        p->xtime=0.0;
        p->printouttime=0.0;

        if(hsmax!=hsmax || nmin<0.0)
        {
            if(p->mpirank==0)
            cout<<endl<<"EMERGENCY STOP  --  SEASTATE: wave height is NaN or the action density negative"<<endl<<endl;

            pprint->print2D(p,e,pgc);
            pgc->final(true);
        }

        // Phase 18: checkpoint and rebalance (nonstationary), the load map (G 45)
        ckpt_step(p,pgc);
    }

    if(p->mpirank==0)
    {
    cout<<endl<<"******************************"<<endl<<endl;

    if(p->plog)
    p->plog->end(p,"finished");

    cout<<"modelled time: "<<p->simtime<<endl;
    cout<<endl;
    }

    // Phase 19: the time of the pipelined sweeps per rank
    if(kba)
    pipeline_report(p,pgc);

    preg->seastate_final(p,e,pgc);

    pgc->final();
}

// Phase 19: the time of the pipelined sweeps of every rank over the run (seastate_implicit::kba_stats), against the time
// of the iterations (transport): on stdout and in REEF3D_SEASTATE_Log/REEF3D_SEASTATE_pipeline.dat
void seastate_f::pipeline_report(lexer *p, ghostcell *pgc)
{
    const seastate_implicit::kba_stats &s = psolv->pipeline_stats();
    const int nv = 10;
    double mine[nv] = {t_trans,double(n_trans),double(s.sweeps),s.wall,s.busy,s.fill,s.idle,s.comm,s.finish,s.exch};
    std::vector<double> all(size_t(nv)*p->mpi_size,0.0);
    MPI_Gather(mine,nv,MPI_DOUBLE,all.data(),nv,MPI_DOUBLE,0,MPI_COMM_WORLD);

    if(p->mpirank!=0)
    return;

    mkdir("./REEF3D_SEASTATE_Log",0777);
    ofstream out("./REEF3D_SEASTATE_Log/REEF3D_SEASTATE_pipeline.dat");
    out<<"# SEASTATE pipelined sweeps (A 781 1), tiles of "<<psolv->pipeline_tile(0)<<" x "<<psolv->pipeline_tile(1)<<" cells; per rank [s] over the run:"<<endl;
    out<<"# the iterations (transport), their number, the sweeps, their wall time, in them the work of the tiles, the waiting before"
       <<" the first tile (the fill), the waiting between the tiles, the messages, the end of the sweeps (the sends completed); the"
       <<" exchanges around the sweeps (halo, fills)"<<endl;
    out<<"# rank \t iterations [s] \t iterations \t sweeps \t sweeps [s] \t work \t fill \t idle \t messages \t end \t exchanges"<<endl;
    out<<setprecision(6);

    const streamsize pr = cout.precision();
    cout<<"SEASTATE pipeline (A 781 1), tiles of "<<psolv->pipeline_tile(0)<<" x "<<psolv->pipeline_tile(1)<<" cells: per rank the share of the iterations"<<endl;
    cout<<"   rank   iterations [s]   sweeps   work   fill   between tiles   messages   end   exchanges"<<endl;
    for(int r=0; r<p->mpi_size; ++r)
    {
        const double *v = &all[size_t(nv)*r];
        out<<r;
        for(int k=0; k<nv; ++k)
        out<<" \t "<<v[k];
        out<<"\n";
        const double T = v[0]>0.0 ? v[0] : 1.0;
        cout<<setprecision(3)<<"   "<<r<<"   "<<v[0]<<"   "<<100.0*v[3]/T<<" %   "<<100.0*v[4]/T<<" %   "<<100.0*v[5]/T<<" %   "
            <<100.0*v[6]/T<<" %   "<<100.0*v[7]/T<<" %   "<<100.0*v[8]/T<<" %   "<<100.0*v[9]/T<<" %"<<endl;
    }
    cout.precision(pr);
}

void seastate_f::step(lexer *p, ghostcell *pgc)
{
    // implicit transport in x, y, sigma, theta with the source terms (if any)
    tnew = p->simtime + dtw;

    if(bser!=nullptr || wmet!=nullptr)
    forcing_update(p,pgc,tnew);

    if(sb!=nullptr)
    surfbeat_boundary(p,tnew);

    const double tt0 = pgc->timer();
    transport(p,pgc);
    t_trans += pgc->timer()-tt0;
    n_trans += iter_done;

    if(sb!=nullptr)
    surfbeat_step(p,pgc);

    parameters(p,pgc);
}

void seastate_f::step_coupled(lexer *p, ghostcell *pgc, double dt)
{
    dtw = dt;

    // the host has written depth, eta, U, V: cells of the initial active set dry and re-wet
    // with the host's water depth; a cell that dries loses its spectrum
    const int nbin = e->grid->nbin;

    IMALOOP
    JMALOOP
    {
    const int w = (e->wet0(i,j)==1 && e->depth(i,j)>=p->A705) ? 1 : 0;

        if(w==0 && e->wet(i,j)==1)
        {
        float *s = e->N->spec(i,j);

        if(s!=nullptr)
        for(int b=0; b<nbin; ++b)
        s[b] = 0.0f;
        }

    e->wet(i,j) = w;
    }

    kinematics(p,pgc,dt);

    // the host calls at the start of its step: the wave field is advanced to the host's time
    tnew = p->simtime;

    if(bser!=nullptr || wmet!=nullptr)
    forcing_update(p,pgc,tnew);

    if(sb!=nullptr)
    surfbeat_boundary(p,tnew);

    transport(p,pgc);

    if(sb!=nullptr)
    surfbeat_step(p,pgc);

    parameters(p,pgc);

    log_step(p);

    if(p->A760>0)
    handover(p,pgc);
}
