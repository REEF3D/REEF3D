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

#include"amr_load.h"
#include"reef3d_mpmd.h"
#include"lexer.h"
#include"reefckpt.h"
#include<mpi.h>
#include<sys/stat.h>
#include<cmath>
#include<iostream>
#include<iomanip>

using namespace std;

double amr_clock::mpi = 0.0;
long   amr_clock::mpi_calls = 0;
double amr_clock::patch = 0.0;
double amr_clock::regrid = 0.0;
double amr_clock::regrid_mpi = 0.0;
int    amr_clock::depth = 0;
int    amr_clock::in_regrid = 0;
double amr_clock::t_patch = 0.0;
double amr_clock::m_patch = 0.0;

void amr_clock::patch_begin()
{
    if(depth++==0 && in_regrid==0)
    {
        t_patch = MPI_Wtime();
        m_patch = mpi;
    }
}

void amr_clock::patch_end()
{
    if(--depth==0 && in_regrid==0)
    patch += (MPI_Wtime()-t_patch)-(mpi-m_patch);
    if(depth<0)
    depth=0;
}

amr_regrid_clock::amr_regrid_clock() : t(MPI_Wtime()), m(amr_clock::mpi)
{
    // a regrid inside a patch scope: the patch time up to here counts, the rest is regrid time
    if(amr_clock::in_regrid==0 && amr_clock::depth>0)
    amr_clock::patch += (t-amr_clock::t_patch)-(amr_clock::mpi-amr_clock::m_patch);
    ++amr_clock::in_regrid;
}

amr_regrid_clock::~amr_regrid_clock()
{
    if(--amr_clock::in_regrid==0)
    {
        const double now = MPI_Wtime();
        amr_clock::regrid += now-t;
        amr_clock::regrid_mpi += amr_clock::mpi-m;
        if(amr_clock::depth>0)
        {
            amr_clock::t_patch = now;
            amr_clock::m_patch = amr_clock::mpi;
        }
    }
}

static double amr_window_last_busy = 0.0;

double amr_window_busy()
{
    return amr_window_last_busy;
}

double amr_window_imbalance()
{
    static bool open = false;
    static double w0 = 0.0, m0 = 0.0;

    const double now = MPI_Wtime();
    const double busy = (now-w0) - (amr_clock::mpi-m0);
    const bool first = !open;
    open = true;

    double mx=0.0, sm=0.0;
    MPI_Allreduce(&busy,&mx,1,MPI_DOUBLE,MPI_MAX,reef3d_comm());
    MPI_Allreduce(&busy,&sm,1,MPI_DOUBLE,MPI_SUM,reef3d_comm());
    int np=1;
    MPI_Comm_size(reef3d_comm(),&np);

    w0 = MPI_Wtime();
    m0 = amr_clock::mpi;
    amr_window_last_busy = first ? 0.0 : busy;

    if(first || sm<=0.0)
    return 0.0;
    return mx/(sm/np);
}

void amr_load::configure(const string &nm)
{
    name = nm;
    tag = nm;
    for(auto &ch : tag)
    if(ch==' ')
    ch='_';
    dir = "./REEF3D_"+tag;
}

void amr_load::step(lexer *p, const function<amr_load_counts()> &counts)
{
    if(p->G43<=0)
    return;

    if(name.empty())
    configure("AMR");

    // the first step opens the window (the initialisation is not counted)
    if(!started)
    {
        const amr_load_counts c = counts();
        started = true;
        steps = 0;
        w0 = MPI_Wtime();
        m0 = amr_clock::mpi;
        p0 = (c.tpatch>=0.0) ? c.tpatch : amr_clock::patch;
        r0 = amr_clock::regrid;
        rm0 = amr_clock::regrid_mpi;
        return;
    }

    ++steps;

    if(p->count%p->G43==0)
    report(p,counts(),steps);
}

void amr_load::report(lexer *p, const amr_load_counts &c, int nsteps)
{
    const int NV = 11;
    const double now = MPI_Wtime();
    const double tp = (c.tpatch>=0.0) ? c.tpatch : amr_clock::patch;

    double d[NV] = {now-w0, amr_clock::mpi-m0, tp-p0, amr_clock::regrid-r0, amr_clock::regrid_mpi-rm0,
                    c.n0, c.n0a, c.np, double(c.npatch), c.wpred, double(nsteps)};

    // the next window starts here (the gather below counts into it)
    w0 = now;
    m0 = amr_clock::mpi;
    p0 = tp;
    r0 = amr_clock::regrid;
    rm0 = amr_clock::regrid_mpi;
    steps = 0;

    const int nr = p->mpi_size;
    const int me = p->mpirank;
    vector<double> all(me==0 ? (size_t)NV*nr : 1);
    MPI_Gather(d,NV,MPI_DOUBLE,&all[0],NV,MPI_DOUBLE,0,reef3d_comm());

    if(me!=0)
    return;

    auto V = [&](int r, int k) { return all[(size_t)r*NV+k]; };

    // busy time: the wall-clock time without the MPI calls
    vector<double> busy(nr);
    double swall=0.0, smpi=0.0, spatch=0.0, sreg=0.0, sbusy=0.0, mbusy=-1.0, spred=0.0, mpred=0.0;
    double scells=0.0, spatches=0.0;
    int slow=0;
    for(int r=0; r<nr; ++r)
    {
        busy[r] = V(r,0)-V(r,1);
        swall += V(r,0);
        smpi += V(r,1);
        spatch += V(r,2);
        sreg += V(r,3);
        sbusy += busy[r];
        if(busy[r]>mbusy)
        {
            mbusy = busy[r];
            slow = r;
        }
        spred += V(r,9);
        mpred = max(mpred,V(r,9));
        scells += V(r,5)+V(r,7);
        spatches += V(r,8);
    }
    const double mwall = swall/nr, mb = sbusy/nr, mp = spred/nr;
    const double L = (mb>0.0) ? mbusy/mb : 1.0;
    const double Lp = (mp>0.0) ? mpred/mp : 1.0;
    const double eff = (mwall>0.0) ? mb/mwall : 0.0;
    const double loss = (mwall>0.0) ? (mbusy-mb)/mwall : 0.0;
    const double fm = (mwall>0.0) ? smpi/nr/mwall : 0.0;
    const double fp = (mwall>0.0) ? spatch/nr/mwall : 0.0;
    const double fr = (mwall>0.0) ? sreg/nr/mwall : 0.0;

    // the phases: patch kernels and level 0 (busy without the patches and the work of the regrids)
    vector<double> lev0(nr);
    double sp=0.0, mxp=0.0, s0=0.0, mx0=0.0;
    for(int r=0; r<nr; ++r)
    {
        lev0[r] = busy[r]-V(r,2)-(V(r,3)-V(r,4));
        sp += V(r,2); mxp = max(mxp,V(r,2));
        s0 += lev0[r]; mx0 = max(mx0,lev0[r]);
    }
    const double Lpatch = (sp>0.0) ? mxp/(sp/nr) : 1.0;
    const double Llev0 = (s0>0.0) ? mx0/(s0/nr) : 1.0;

    // the measured cost of a patch cell in level-0 cells: busy = a n0 + b np over the ranks
    double wfit = -1.0;
    {
        double S00=0.0, S01=0.0, S11=0.0, B0=0.0, B1=0.0;
        for(int r=0; r<nr; ++r)
        {
            const double n0 = V(r,5), np = V(r,7);
            S00 += n0*n0; S01 += n0*np; S11 += np*np;
            B0 += n0*busy[r]; B1 += np*busy[r];
        }
        const double det = S00*S11-S01*S01;
        if(nr>=3 && S11>0.0 && det>1.0e-6*S00*S11)
        {
            const double a = (B0*S11-B1*S01)/det;
            const double b = (S00*B1-S01*B0)/det;
            if(a>0.0)
            wfit = b/a;
        }
    }

    // correlation of the predicted work and the busy time over the ranks
    double rcorr = 0.0;
    {
        double sxy=0.0, sxx=0.0, syy=0.0;
        for(int r=0; r<nr; ++r)
        {
            const double x = V(r,9)-mp, y = busy[r]-mb;
            sxy += x*y; sxx += x*x; syy += y*y;
        }
        if(sxx>0.0 && syy>0.0)
        rcorr = sxy/sqrt(sxx*syy);
    }

    if(!sout.is_open())
    {
        mkdir(dir.c_str(),0777);
        // I 45: a restart continues the files
        const bool h1 = reefckpt_open(sout,(dir+"/REEF3D_"+tag+"_load.dat").c_str(),p->I45);
        const bool h2 = reefckpt_open(rout,(dir+"/REEF3D_"+tag+"_load_ranks.dat").c_str(),p->I45);
        if(h1)
        {
        sout<<"# "<<name<<" load statistics (G 43 "<<p->G43<<"), "<<nr<<" ranks; busy = wall - MPI time per rank"<<endl;
        sout<<"# count \t simtime \t steps \t wall per step [s] \t busy max/mean \t predicted max/mean \t parallel efficiency mean(busy)/wall"
              " \t imbalance loss (max-mean busy)/wall \t MPI share \t patch share \t regrid share \t patch cell cost (fit, level-0 cells)"
              " \t r(predicted,busy) \t slowest rank \t patches \t cells \t patch time max/mean \t level-0 time max/mean"<<endl;
        }
        if(h2)
        {
        rout<<"# "<<name<<" load per rank (G 43 "<<p->G43<<"); times in s over the steps since the last report"<<endl;
        rout<<"# count \t rank \t level-0 cells \t not solid \t patch cells \t patches \t predicted work \t wall \t MPI \t patch \t regrid \t busy \t busy/mean \t level 0 (busy - patch - regrid work)"<<endl;
        }
    }

    sout<<p->count<<" \t "<<setprecision(10)<<p->simtime<<" \t "<<nsteps<<" \t "<<setprecision(5)<<mwall/max(nsteps,1)
        <<" \t "<<L<<" \t "<<Lp<<" \t "<<eff<<" \t "<<loss<<" \t "<<fm<<" \t "<<fp<<" \t "<<fr<<" \t "<<wfit<<" \t "<<rcorr
        <<" \t "<<slow<<" \t "<<long(spatches)<<" \t "<<setprecision(12)<<scells<<" \t "<<setprecision(5)<<Lpatch<<" \t "<<Llev0<<endl;

    for(int r=0; r<nr; ++r)
    rout<<p->count<<" \t "<<r<<" \t "<<setprecision(10)<<V(r,5)<<" \t "<<V(r,6)<<" \t "<<V(r,7)<<" \t "<<long(V(r,8))<<" \t "<<V(r,9)
        <<" \t "<<setprecision(5)<<V(r,0)<<" \t "<<V(r,1)<<" \t "<<V(r,2)<<" \t "<<V(r,3)<<" \t "<<busy[r]<<" \t "<<(mb>0.0 ? busy[r]/mb : 1.0)<<" \t "<<lev0[r]<<endl;

    sout.flush();
    rout.flush();

    const streamsize pr = cout.precision();
    cout<<name<<" load over "<<nsteps<<" steps on "<<nr<<" ranks: busy max/mean "<<setprecision(3)<<L<<" (predicted "<<Lp
        <<"), by phase: patches "<<Lpatch<<", level 0 "<<Llev0<<"; efficiency "<<eff<<", lost to imbalance "<<loss<<" of the wall time; MPI "<<fm
        <<", patches "<<fp<<", regrid "<<fr<<"; slowest rank "<<slow;
    if(wfit>0.0)
    cout<<"; a patch cell costs "<<wfit<<" level-0 cells";
    if(nr>=3)
    cout<<"; r(predicted,busy) "<<rcorr;
    cout<<endl;
    cout.precision(pr);
}
