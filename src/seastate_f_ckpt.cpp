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

//  Phase 18: SEASTATE checkpoints on any decomposition (P 48, I 45) and the rebalance on the measured
//  work of the cells (G 45, G 46), for the run part of reefrun_ckpt.h.
//
//  The state: the spectra of the active cells of all grids (level 0 and the patches of the mesh refinement;
//  the dry cells in the blocks of the block-sparse storage hold no energy, and the blocks depend on the ranks),
//  with the active ranges of the spectral sparsity (A 795) and, within a stationary step, the last
//  change of Hs of the convergence test (A 799 1), as records by the global index of the cell on its
//  level (rank blob ss_cells: level, I, J, the change of Hs, the ranges, the spectrum as floats, two
//  per double, bitwise).  A restart takes the records of the pieces whose level-0 box meets its own
//  (the patches are cut at the rank boxes), the restriction, the level-0 halo and the fills follow
//  (seastate_amr::restored).  A stationary step continues at the iteration of the checkpoint (the
//  checkpoints and the checks of the load are counted in iterations there).  On the same
//  decomposition the restart continues bitwise; on another one as well, since the sweeps give the
//  result of one rank (pipelined, A 781 1) -- the lagged sweeps differ as they differ between ranks.
//  Only the latest checkpoint of a run is kept (a spectrum per cell is large).
//
//  The work: the time of the solves of every cell (seastate_implicit::work_field), per level-0 column
//  of the rank with the patch cells over it.  The load of a rank is the sum over its cells (the waits
//  of the pipeline are not in it: they spin outside MPI), the load map the work per iteration of every
//  level-0 column [microseconds]: REEF3D_SEASTATE_Log/REEF3D_SEASTATE_loadmap.dat, which DIVEMesh
//  (M 20 5, M 35 1) and the rebalance restart (M 41 1) take.
//
//  Not covered: the surfbeat model (A 770 1), runs coupled to a host model (A 750).

#include"seastate_f.h"
#include"reef3d_mpmd.h"
#include"fdm_seastate.h"
#include"seastate_grid.h"
#include"seastate_store.h"
#include"seastate_exchange.h"
#include"seastate_implicit.h"
#include"seastate_amr.h"
#include"slice4.h"
#include"sliceint4.h"
#include"lexer.h"
#include"ghostcell.h"
#include"reefckpt.h"
#include"reefrun_ckpt.h"
#include<algorithm>
#include<cmath>
#include<cstring>
#include<iomanip>
#include<iostream>
#include<sstream>
#include<unordered_map>
#include<mpi.h>
#include<sys/stat.h>

using namespace std;

namespace
{
// a cell of level lev at the global index (I, J) of its level
inline long long ckey(int lev, int I, int J)
{
    return ((long long)lev<<56) | ((long long)(I+(1<<27))<<28) | (long long)(J+(1<<27));
}
}

string seastate_f::ckpt_unsupported(lexer *p) const
{
    if(coupled)
    return "SEASTATE coupled to a host model (A 750) is not in the checkpoints";
    if(p->A770==1)
    return "the SEASTATE surfbeat model (A 770 1) is not in the checkpoints yet";
    return "";
}

reefrun_ckpt& seastate_f::ckpt(lexer *p, ghostcell *pgc)
{
    if(rck!=nullptr)
    return *rck;

    reefrun_ckpt::hooks h;

    h.unsupported = [this,p](){ return ckpt_unsupported(p); };
    h.write = [this,p](reefckpt &ck){ ckpt_write(p,ck); };
    h.read = [this,p,pgc](reefckpt &ck, string &msg){ return ckpt_read(p,pgc,ck,msg); };

    // the outputs of the patches (the VTP of level 0 counts with p->printcount)
    h.printer_state = [this](vector<double> &v){ v.push_back(pamr!=nullptr ? double(pamr->print_count()) : 0.0); };
    h.printer_restore = [this](const vector<double> &v){ if(pamr!=nullptr && !v.empty()) pamr->print_count(int(v[0])); };

    h.flush = [this]()
    {
        integral.flush();
        convlog.flush();
        if(pamr!=nullptr)
        pamr->flush_logs();
    };

    // the load map of the window the rebalance measured (G 46), and the load of the ranks in it
    h.loadmap = [this](lexer *pp){ loadmap_write(pp,wwin,nwin,"the last G 47 window"); };
    h.imbalance = [this,p,pgc](){ return work_imbalance(p,pgc); };

    // the load max/mean of a partition by the load map of the window (rank 0)
    h.predict = [this,p](const vector<int> &nb)
    {
        const int GNY = p->gknoy, np = int(nb.size()/6);
        if(lmap0.size()!=size_t(p->gknox)*GNY || np<1)
        return 0.0;
        double mx = 0.0, sum = 0.0;
        for(int r=0; r<np; ++r)
        {
            double s = 0.0;
            for(int a=0; a<nb[6*r+3]; ++a)
            for(int b=0; b<nb[6*r+4]; ++b)
            s += lmap0[size_t(nb[6*r]+a)*GNY+(nb[6*r+1]+b)];
            mx = std::max(mx,s);
            sum += s;
        }
        return sum>0.0 ? mx*double(np)/sum : 0.0;
    };

    // the gain of the lower load over the iterations left (stationary: estimated, nonstationary: the steps left) against
    // the restart: the set-up of this run, the checkpoint written and read (at the speed of the last one, else 100 MB/s)
    h.pays = [this,p](double L, double Lnew, string &why)
    {
        const double gain_it = (nwin>0.0 && L>0.0) ? busy_max/nwin*(1.0-Lnew/L) : 0.0;
        double rem;
        bool est = false;
        if(p->A700==2)
        {
            // stationary: the iterations left by the trend of the last G 47 iterations (since the start of the step or the
            // restart): the largest change of Hs (geometric) to A 708, the percentage of the cells within A 708 (linear) to
            // A 709; without a trend as many as so far; at most up to A 707
            rem = double(std::max(iter_max-k_now,0));
            const double lindy = std::min(rem,double(k_now));
            const size_t n = hconv.size(), w = std::min(size_t(std::max(p->G47,1)),n>0 ? n-1 : 0);
            if(w>0)
            {
                const double c1 = hconv[n-1], c0 = hconv[n-1-w];
                if(c1>0.0 && c0>c1 && p->A708>0.0)
                {
                    rem = std::min(rem,std::max(0.0,std::log(p->A708/c1)/std::log(std::pow(c1/c0,1.0/double(w)))));
                    est = true;
                }
                const double q1 = hperc[n-1], q0 = hperc[n-1-w];
                if(p->A709<100.0 && q1>q0)
                {
                    rem = std::min(rem,std::max(0.0,(p->A709-q1)*double(w)/(q1-q0)));
                    est = true;
                }
            }
            // no trend yet (the change of Hs and the percentage not converging): as many iterations as so far
            if(!est)
            rem = lindy;
        }
        else
        {
            double st = double(p->N45-p->count);
            if(p->dt>0.0)
            st = std::min(st,(p->N41-p->simtime)/p->dt);
            rem = std::max(st,0.0)*double(iter_max);
        }
        const double bw = (rck!=nullptr && rck->write_s>0.0) ? state_bytes/rck->write_s : 1.0e8;
        const double tck = bw>0.0 ? state_bytes/bw : 0.0;
        const double cost = t_ini+2.0*tck, gain = gain_it*rem;
        ostringstream o;
        o<<setprecision(3)<<"the restart costs about "<<cost<<" s (the set-up "<<t_ini<<" s, the checkpoint of "<<state_bytes/1.0e9
         <<" GB written and read "<<2.0*tck<<" s), the lower load saves about "<<gain_it<<" s per iteration over "
         <<"about "<<long(std::ceil(rem))<<(long(std::ceil(rem))==1 ? " iteration" : " iterations")<<(p->A700!=2 ? " (the steps left)" : (est ? " (estimated from the convergence)" : " (as many as so far)"))<<" ("<<gain<<" s)";
        why = o.str();
        return gain>cost;
    };

    h.keep_last = true;

    rck = new reefrun_ckpt(p,pgc,std::move(h));
    return *rck;
}

// ------------------------------------------------------------------------------------------- the checks

void seastate_f::ckpt_iteration(lexer *p, int k)
{
    if(rck==nullptr || (p->P48<=0 && p->G46<=0.0))
    return;

    // not after the last iteration of the step: the step ends with it
    rck->sub = k;
    k_now = k;
    if(p->P48>0 && k%p->P48==0 && k<iter_max)
    rck->write();
    // a stop signal (SIGTERM, SIGUSR1): within the step too (a stationary step can be long)
    if(p->P48>0 && k<iter_max)
    rck->signal_stop();
    if(p->G46>0.0)
    rck->rebalance_check_iter(k,iter_max);
    rck->sub = 0;
}

void seastate_f::ckpt_step(lexer *p, ghostcell *pgc)
{
    // G 45: the load map of the run so far every G 45 steps
    if(wk0!=nullptr && p->G45>0 && p->count%p->G45==0)
    {
        vector<double> col;
        work_columns(p,col,false);
        if(wrun.size()==col.size())
        for(size_t n=0; n<col.size(); ++n)
        col[n] += wrun[n];
        loadmap_write(p,col,nrun+double(wit),"the run");
    }

    // stationary: the checkpoints and the checks of the load within the step (ckpt_iteration)
    if(rck==nullptr || p->A700==2)
    return;

    rck->step();
}

// ------------------------------------------------------------------------------------------- the work

void seastate_f::work_columns(lexer *p, vector<double> &col, bool reset)
{
    col.assign(size_t(p->knox)*p->knoy,0.0);

    if(wk0!=nullptr)
    for(int ii=0; ii<p->knox; ++ii)
    for(int jj=0; jj<p->knoy; ++jj)
    {
        col[size_t(ii)*p->knoy+jj] += (*wk0)(ii,jj);
        if(reset)
        (*wk0)(ii,jj) = 0.0;
    }

    if(pamr!=nullptr)
    pamr->work_columns(p,col,reset);
}

// the load max/mean of the ranks in the window since the last call (collective); the window is kept for the load map
double seastate_f::work_imbalance(lexer *p, ghostcell *pgc)
{
    vector<double> col;
    work_columns(p,col,true);

    if(wrun.size()!=col.size())
    wrun.assign(col.size(),0.0);
    double busy = 0.0;
    for(size_t n=0; n<col.size(); ++n)
    {
        wrun[n] += col[n];
        busy += col[n];
    }
    wwin.swap(col);
    nwin = double(wit);
    nrun += double(wit);
    wit = 0;

    const double mx = pgc->globalmax(busy);
    const double mean = pgc->globalsum(busy)/double(p->mpi_size);
    busy_max = mx;
    busy_mean = mean;

    // the bytes of the state in a checkpoint (the active cells of all grids; the grids are static: once)
    if(state_bytes<=0.0)
    {
        double nc = 0.0;
        for(int ii=0; ii<p->knox; ++ii)
        for(int jj=0; jj<p->knoy; ++jj)
        if(e->N->spec(ii,jj)!=nullptr && e->wet(ii,jj)==1)
        nc += 1.0;
        if(pamr!=nullptr)
        for(int id=0; id<pamr->npatch(); ++id)
        {
            const seastate_amr::gview G = pamr->patch_view(id);
            for(int ii=G.ox; ii<G.ox+G.nx; ++ii)
            for(int jj=G.oy; jj<G.oy+G.ny; ++jj)
            if(G.e->N->spec(ii,jj)!=nullptr && G.e->wet(ii,jj)==1)
            nc += 1.0;
        }
        state_bytes = pgc->globalsum(nc)*(8.0*double(4+psolv->ranges_size()+(e->grid->nbin+1)/2));
    }

    if(p->mpirank==0 && mean>0.0)
    {
        const streamsize pr = cout.precision();
        cout<<"SEASTATE work of the cells in the window: "<<setprecision(4)<<mean*double(p->mpi_size)<<" s on all ranks, the heaviest rank "
            <<mx<<" s ("<<long(nwin)<<" iterations)"<<endl;
        cout.precision(pr);
    }

    return mean>0.0 ? mx/mean : 0.0;
}

// the load map (collective): per level-0 column of the domain the work per iteration [microseconds]
void seastate_f::loadmap_write(lexer *p, const vector<double> &col, double iterations, const char *what)
{
    const int np = p->mpi_size;
    int mybox[4] = {p->origin_i,p->origin_j,p->knox,p->knoy};
    vector<int> boxes(4*size_t(np));
    MPI_Gather(mybox,4,MPI_INT,boxes.data(),4,MPI_INT,0,reef3d_comm());

    vector<double> mine(size_t(p->knox)*p->knoy,0.0);
    const double f = iterations>0.0 ? 1.0e6/iterations : 0.0;
    if(col.size()==mine.size())
    for(size_t n=0; n<mine.size(); ++n)
    mine[n] = f*col[n];

    vector<int> cnt(np), off(np,0);
    if(p->mpirank==0)
    {
        for(int r=0; r<np; ++r)
        cnt[r] = boxes[4*r+2]*boxes[4*r+3];
        for(int r=1; r<np; ++r)
        off[r] = off[r-1]+cnt[r-1];
    }
    vector<double> all(p->mpirank==0 ? std::max(size_t(off[np-1]+cnt[np-1]),size_t(1)) : 1);
    MPI_Gatherv(mine.data(),(int)mine.size(),MPI_DOUBLE,all.data(),cnt.data(),off.data(),MPI_DOUBLE,0,reef3d_comm());

    if(p->mpirank!=0)
    return;

    const int GNX = p->gknox, GNY = p->gknoy;
    vector<double> map(size_t(GNX)*GNY,0.0);
    for(int r=0; r<np; ++r)
    {
        const int oi = boxes[4*r], oj = boxes[4*r+1], ni = boxes[4*r+2], nj = boxes[4*r+3];
        for(int a=0; a<ni; ++a)
        for(int b=0; b<nj; ++b)
        map[size_t(oi+a)*GNY+(oj+b)] = all[off[r]+size_t(a)*nj+b];
    }

    lmap0 = map;

    mkdir("./REEF3D_SEASTATE_Log",0777);
    ofstream out("./REEF3D_SEASTATE_Log/REEF3D_SEASTATE_loadmap.dat");
    out<<"# SEASTATE load map (G 45, G 46): the measured work of the cells per iteration of every level-0 column [microseconds], "
       <<"the patch cells over it included, over "<<what<<" ("<<long(iterations)<<" iterations)"<<endl;
    out<<"# for the weighted partition of DIVEMesh: M 20 5, M 35 1 (copy as loadmap.dat or run DIVEMesh in this case folder)"<<endl;
    out<<GNX<<" "<<GNY<<endl;
    out<<setprecision(6);
    for(int I=0; I<GNX; ++I)
    {
        for(int J=0; J<GNY; ++J)
        out<<map[size_t(I)*GNY+J]<<(J+1<GNY ? " " : "");
        out<<"\n";
    }
}

// ------------------------------------------------------------------------------------------- checkpoints

void seastate_f::ckpt_write(lexer *p, reefckpt &ck)
{
    const int nbin = e->grid->nbin;
    const int nr = psolv->ranges_size();
    const int nl = psolv->lag_size();
    const int ns = (nbin+1)/2;
    const int k = (rck!=nullptr) ? rck->sub : 0;
    const bool mid = (k>0 && p->A700==2);

    ck.scalar("ss_iter",double(mid ? k : 0));
    ck.scalar("ss_nbin",double(nbin));
    ck.scalar("ss_nr",double(nr));
    ck.scalar("ss_nl",double(nl));     // Phase 23: the state of the wind input per cell (A 768 1)

    // the last change of Hs of the leaf cells (within a stationary step): level 0 by cell, the patches by cell
    vector<double> d0(size_t(p->knox)*p->knoy,0.0);
    vector<vector<double>> dp(pamr!=nullptr ? pamr->npatch() : 0);
    if(mid)
    {
        size_t n = 0;
        if(pamr==nullptr)
        {
            if(dprev.size()==d0.size())
            d0 = dprev;
        }
        else
        {
            SLICELOOP4
            if(e->wet(i,j)==1 && pamr->level0_covered(i,j)==0)
            {
            if(n<dprev.size())
            d0[size_t(i)*p->knoy+j] = dprev[n];
            ++n;
            }

            for(int id=0; id<pamr->npatch(); ++id)
            {
            const seastate_amr::gview G = pamr->patch_view(id);
            lexer *pp = G.pp;
            dp[id].assign(size_t(pp->imax)*pp->jmax,0.0);
                for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
                for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
                if(G.e->wet(ii,jj)==1 && ii>=G.ox && ii<G.ox+G.nx && jj>=G.ox && jj<G.ox+G.ny && (*G.cov)(ii,jj)==0)
                {
                if(n<dprev.size())
                dp[id][size_t(ii-pp->imin)*pp->jmax+(jj-pp->jmin)] = dprev[n];
                ++n;
                }
            }
        }
    }

    vector<double> b;
    vector<double> r(std::max(nr,1)), lg(std::max(nl,1));
    auto put = [&](int lev, int I, int J, double d, const float *N, seastate_implicit *s, lexer *pp, int ii, int jj)
    {
        b.push_back(double(lev));
        b.push_back(double(I));
        b.push_back(double(J));
        b.push_back(d);
        if(nr>0)
        {
            s->ranges_get(pp,ii,jj,r.data());
            b.insert(b.end(),r.begin(),r.begin()+nr);
        }
        if(nl>0)
        {
            s->lag_get(pp,ii,jj,lg.data());
            b.insert(b.end(),lg.begin(),lg.begin()+nl);
        }
        const size_t at = b.size();
        b.resize(at+ns,0.0);
        memcpy(&b[at],N,sizeof(float)*size_t(nbin));
    };

    for(int ii=0; ii<p->knox; ++ii)
    for(int jj=0; jj<p->knoy; ++jj)
    {
        const float *N = e->N->spec(ii,jj);
        if(N!=nullptr && e->wet(ii,jj)==1)
        put(0,ii+p->origin_i,jj+p->origin_j,d0[size_t(ii)*p->knoy+jj],N,psolv,p,ii,jj);
    }

    if(pamr!=nullptr)
    for(int id=0; id<pamr->npatch(); ++id)
    {
        const seastate_amr::gview G = pamr->patch_view(id);
        lexer *pp = G.pp;
        for(int ii=G.ox; ii<G.ox+G.nx; ++ii)
        for(int jj=G.oy; jj<G.oy+G.ny; ++jj)
        {
            const float *N = G.e->N->spec(ii,jj);
            if(N!=nullptr && G.e->wet(ii,jj)==1)
            put(G.lev,ii-G.ox+G.I0,jj-G.oy+G.J0,mid ? dp[id][size_t(ii-pp->imin)*pp->jmax+(jj-pp->jmin)] : 0.0,N,G.solv,pp,ii,jj);
        }
    }

    ck.rank_blob("ss_cells",b);
}

bool seastate_f::ckpt_read(lexer *p, ghostcell *pgc, reefckpt &ck, string &msg)
{
    const int nbin = e->grid->nbin;
    const int nr = psolv->ranges_size();
    const int nl = psolv->lag_size();
    const int ns = (nbin+1)/2;

    if(!ck.has("ss_nbin"))
    {
        msg = "the checkpoint holds no SEASTATE state";
        return false;
    }
    if(int(ck.get_scalar("ss_nbin"))!=nbin)
    {
        msg = "the spectral grid of the checkpoint ("+to_string(int(ck.get_scalar("ss_nbin")))+" bins) differs from A 701 - A 703 ("+to_string(nbin)+")";
        return false;
    }
    const int nrc = int(ck.get_scalar("ss_nr"));
    const int nlc = ck.has("ss_nl") ? int(ck.get_scalar("ss_nl")) : 0;
    const int rec = 4+nrc+nlc+ns;
    // the state of the wind input (A 768 1): only if both runs keep it, else it starts from zero
    const bool lags = (nl>0 && nlc==nl);
    // the ranges: only with the same sparsity, else they start again (as at the first visit of a cell)
    const bool ranges = (nr>0 && nrc==nr);

    // the cells of this rank: grid (-1 level 0, else the patch), local index; restored or not
    struct dst {int g, ii, jj; bool got;};
    unordered_map<long long,dst> mine;
    mine.reserve(size_t(p->knox)*p->knoy*2);
    long ncells = 0;

    for(int ii=0; ii<p->knox; ++ii)
    for(int jj=0; jj<p->knoy; ++jj)
    if(e->N->spec(ii,jj)!=nullptr && e->wet(ii,jj)==1)
    {
        mine[ckey(0,ii+p->origin_i,jj+p->origin_j)] = dst{-1,ii,jj,false};
        ++ncells;
    }

    int maxl = 0;
    if(pamr!=nullptr)
    for(int id=0; id<pamr->npatch(); ++id)
    {
        const seastate_amr::gview G = pamr->patch_view(id);
        maxl = std::max(maxl,G.lev);
        for(int ii=G.ox; ii<G.ox+G.nx; ++ii)
        for(int jj=G.oy; jj<G.oy+G.ny; ++jj)
        if(G.e->N->spec(ii,jj)!=nullptr && G.e->wet(ii,jj)==1)
        {
            mine[ckey(G.lev,ii-G.ox+G.I0,jj-G.oy+G.J0)] = dst{id,ii,jj,false};
            ++ncells;
        }
    }

    auto grid = [&](const dst &d, lexer *&pp, fdm_seastate *&ee, seastate_implicit *&s)
    {
        if(d.g<0)
        {
            pp = p; ee = e; s = psolv;
            return;
        }
        const seastate_amr::gview G = pamr->patch_view(d.g);
        pp = G.pp; ee = G.e; s = G.solv;
    };

    // the last change of Hs per cell (within a stationary step)
    vector<double> d0(size_t(p->knox)*p->knoy,0.0);
    vector<vector<double>> dp(pamr!=nullptr ? pamr->npatch() : 0);
    for(int id=0; id<(int)dp.size(); ++id)
    {
        lexer *pp = pamr->patch_view(id).pp;
        dp[id].assign(size_t(pp->imax)*pp->jmax,0.0);
    }

    // the pieces whose level-0 box meets the box of this rank (all of them if patches may lie on other ranks, G 40)
    vector<vector<double>> blobs;
    vector<pair<dst*,const double*>> found;     // the records taken (for the spectra once more after the restriction)
    for(int r=0; r<ck.pieces(); ++r)
    {
        int oi, oj, ni, nj;
        ck.piece_box(r,oi,oj,ni,nj);
        const bool meets = oi<p->origin_i+p->knox && oi+ni>p->origin_i && oj<p->origin_j+p->knoy && oj+nj>p->origin_j;
        if(!meets && p->G40<=0)
        continue;

        blobs.push_back(ck.get_rank_blob("ss_cells",r));
        const vector<double> &bl = blobs.back();
        for(size_t q=0; q+rec<=bl.size(); q+=rec)
        {
            auto it = mine.find(ckey(int(bl[q]),int(bl[q+1]),int(bl[q+2])));
            if(it==mine.end())
            continue;
            dst &d = it->second;
            lexer *pp; fdm_seastate *ee; seastate_implicit *s;
            grid(d,pp,ee,s);
            memcpy(ee->N->spec(d.ii,d.jj),&bl[q+4+nrc+nlc],sizeof(float)*size_t(nbin));
            found.push_back({&d,&bl[q+4+nrc+nlc]});
            if(ranges)
            s->ranges_put(pp,d.ii,d.jj,&bl[q+4]);
            if(lags)
            s->lag_put(pp,d.ii,d.jj,&bl[q+4+nrc]);
            if(d.g<0)
            d0[size_t(d.ii)*p->knoy+d.jj] = bl[q+3];
            else
            dp[d.g][size_t(d.ii-pp->imin)*pp->jmax+(d.jj-pp->jmin)] = bl[q+3];
            d.got = true;
        }
    }

    // level 0 must be complete; a patch cell missing (a refinement cut differently) takes the spectrum of its parent
    long miss0 = 0, missp = 0;
    for(auto &m : mine)
    if(!m.second.got && m.second.g<0)
    ++miss0;
    miss0 = long(pgc->globalsum(double(miss0)));
    if(miss0>0)
    {
        msg = to_string(miss0)+" level-0 cells of SEASTATE are not in the checkpoint";
        return false;
    }
    for(int l=1; l<=maxl; ++l)
    for(auto &m : mine)
    {
        dst &d = m.second;
        if(d.got || d.g<0 || pamr->patch_view(d.g).lev!=l)
        continue;
        const seastate_amr::gview G = pamr->patch_view(d.g);
        const int I = d.ii-G.ox+G.I0, J = d.jj-G.oy+G.J0;
        auto it = mine.find(ckey(l-1,I>>1,J>>1));
        if(it!=mine.end() && it->second.got)
        {
            lexer *pp; fdm_seastate *ee; seastate_implicit *s;
            grid(it->second,pp,ee,s);
            memcpy(G.e->N->spec(d.ii,d.jj),ee->N->spec(it->second.ii,it->second.jj),sizeof(float)*size_t(nbin));
            d.got = true;
        }
        ++missp;
    }
    missp = long(pgc->globalsum(double(missp)));

    // the faces of the covered cells, the halo and the ring cells from the spectra of the interiors; the covered cells keep
    // the spectra of the checkpoint (the restriction of a FAS correction rounds otherwise than restrict_all)
    if(pamr!=nullptr)
    pamr->restored(p,pgc,[&]()
    {
        for(const auto &fd : found)
        {
            lexer *pp; fdm_seastate *ee; seastate_implicit *s;
            grid(*fd.first,pp,ee,s);
            memcpy(ee->N->spec(fd.first->ii,fd.first->jj),fd.second,sizeof(float)*size_t(nbin));
        }
    });
    else
    pex->start(p,pgc,*e->N);

    // within a stationary step: the iteration, the last changes of Hs of the leaf cells in the order of the convergence test
    ss_resume = (p->A700==2) ? int(ck.get_scalar("ss_iter",0.0)) : 0;
    dprev.clear();
    if(ss_resume>0)
    {
        if(pamr==nullptr)
        dprev = d0;
        else
        {
            SLICELOOP4
            if(e->wet(i,j)==1 && pamr->level0_covered(i,j)==0)
            dprev.push_back(d0[size_t(i)*p->knoy+j]);

            for(int id=0; id<pamr->npatch(); ++id)
            {
            const seastate_amr::gview G = pamr->patch_view(id);
            lexer *pp = G.pp;
                for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
                for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
                if(G.e->wet(ii,jj)==1 && ii>=G.ox && ii<G.ox+G.nx && jj>=G.ox && jj<G.ox+G.ny && (*G.cov)(ii,jj)==0)
                dprev.push_back(dp[id][size_t(ii-pp->imin)*pp->jmax+(jj-pp->jmin)]);
            }
        }
    }

    const double nall = pgc->globalsum(double(ncells));
    if(p->mpirank==0)
    {
        cout<<"SEASTATE restart: the spectra of "<<long(nall)<<" cells"<<(pamr!=nullptr ? " of all grids" : "")
            <<(ranges ? " with the active ranges of the sparsity" : "");
        if(missp>0)
        cout<<", "<<missp<<" patch cells not in the checkpoint from their parents";
        if(ss_resume>0)
        cout<<"; the stationary step "<<int(ck.get_scalar("count"))<<" continues after iteration "<<ss_resume;
        cout<<endl;
    }

    return true;
}
