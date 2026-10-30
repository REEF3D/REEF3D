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

#include"seastate_implicit.h"
#include"seastate_exchange.h"
#include"seastate_grid.h"
#include"seastate_store.h"
#include"seastate_dispersion.h"
#include"seastate_source.h"
#include"seastate_obstacle.h"
#include"fdm_seastate.h"
#include"sliceint.h"
#include"slice.h"
#include"lexer.h"
#include"ghostcell.h"
#include<algorithm>
#include<atomic>
#include<chrono>
#include<cstdio>
#include<cstdlib>
#include<condition_variable>
#include<memory>
#include<mutex>
#include<thread>
#include<cmath>
#include<mpi.h>

namespace
{
using kclk = std::chrono::steady_clock;

// Phase 12: persistent worker threads for the threaded sweeps (A 798). A sweep hands a job to the
// workers and works its own share; the workers stay awake for a short time after a job (the next
// sweep of a patch follows within microseconds), then sleep until the next one.
class seastate_pool
{
public:
    static seastate_pool &get()
    {
        static seastate_pool p;
        return p;
    }

    // f(t) for t = 0 .. nt-1: t = 0 on the calling thread, the others on the workers
    void run(int nt, const std::function<void(int)> &f)
    {
        if(nt<=1)
        {
        f(0);
        return;
        }

        ensure(nt-1);

        {
        std::lock_guard<std::mutex> lk(mx);
        job = &f;
        njob = nt;
        pending.store(nt-1,std::memory_order_relaxed);
        gen.fetch_add(1,std::memory_order_release);
        }
        cv.notify_all();

        f(0);

        for(int k=0; pending.load(std::memory_order_acquire)>0; ++k)
        if(k>2000)
        std::this_thread::yield();
    }

    ~seastate_pool()
    {
        {
        std::lock_guard<std::mutex> lk(mx);
        quit = true;
        gen.fetch_add(1,std::memory_order_release);
        }
        cv.notify_all();
        for(auto &t : th)
        t.join();
    }

private:
    void ensure(int n)
    {
        // the generation before the job that needs the new workers: they take that job too
        const long g0 = gen.load(std::memory_order_acquire);
        while(int(th.size())<n)
        {
        const int t = int(th.size())+1;
        th.emplace_back([this,t,g0]() {loop(t,g0);});
        }
    }

    void loop(int t, long seen)
    {
        for(;;)
        {
            // short spin for the next job, then sleep
            const auto t0 = std::chrono::steady_clock::now();
            while(gen.load(std::memory_order_acquire)==seen && std::chrono::steady_clock::now()-t0<std::chrono::microseconds(200))
            std::this_thread::yield();

            if(gen.load(std::memory_order_acquire)==seen)
            {
            std::unique_lock<std::mutex> lk(mx);
            cv.wait(lk,[&]() {return gen.load(std::memory_order_acquire)!=seen;});
            }

            const std::function<void(int)> *f;
            int n;
            {
            std::lock_guard<std::mutex> lk(mx);
            if(quit)
            return;
            seen = gen.load(std::memory_order_acquire);
            f = job;
            n = njob;
            }

            if(t<n)
            {
            (*f)(t);
            pending.fetch_sub(1,std::memory_order_acq_rel);
            }
        }
    }

    std::vector<std::thread> th;
    std::mutex mx;
    std::condition_variable cv;
    const std::function<void(int)> *job = nullptr;
    int njob = 0;
    std::atomic<long> gen{0};
    std::atomic<int> pending{0};
    bool quit = false;
};
}

seastate_implicit::seastate_implicit(lexer *p, fdm_seastate *e) : src(nullptr), wU(nullptr), wD(nullptr), second(false)
{
    for(int k=0; k<4; ++k)
    Nbs[k] = Nbs0[k] = nullptr;

    const seastate_grid &g = *e->grid;

    nsig = g.nsig;
    ndir = g.ndir;

    P.assign(g.nbin,0.0);
    D.assign(g.nbin,0.0);

    m0.assign(4,ndir);
    m1.assign(4,-1);

    // the directions of each quadrant (with a computational sector, A 717: the active ones)
    for(int m=0; m<ndir; ++m)
    if(g.act.empty() || g.act[m])
    {
    m0[g.quad[m]] = std::min(m0[g.quad[m]],m);
    m1[g.quad[m]] = std::max(m1[g.quad[m]],m);
    }

    Ar.assign(nsig,0.0);
    qc.assign(ndir,0.0);
    rth.assign(ndir,0.0);
    cdp.assign(size_t(nsig)*ndir,0.0);
    cdm.assign(size_t(nsig)*ndir,0.0);
    qs.assign(ndir,0.0);
    cur.assign(ndir,0.0);
    tAp.assign(ndir,0.0);
    tCp.assign(ndir,0.0);
    tAm.assign(ndir,0.0);
    tCm.assign(ndir,0.0);
    csg.assign(size_t(nsig)*ndir,0.0);
    zero.assign(g.nbin,0.0f);
    zerod.assign(g.nbin,0.0);
    for(vector<double> *v : {&tdi,&tla,&tup,&trh,&tsi,&trd,&tcp,&tdp,&tso})
    v->assign(g.nbin,0.0);
    tok.assign(nsig,1);
    rsum.assign(nsig,0.0);
    wlo.assign(nsig,0);
    whi.assign(nsig,-1);
    thr.assign(nsig,0.0);
    la.assign(ndir,0.0);
    di.assign(ndir,0.0);
    up.assign(ndir,0.0);
    rhs.assign(ndir,0.0);
    sol.assign(ndir,0.0);
    cp.assign(ndir,0.0);
    dp.assign(ndir,0.0);
}

void seastate_implicit::iterate(lexer *p, ghostcell *pgc, fdm_seastate *e, seastate_exchange *pex, const seastate_store *N0,
                                double rdt, const vector<float> &Nb, const int side[4], bool refraction, bool fshift)
{
    kex = pex;

    for(int q=0; q<4; ++q)
    {
    sweep(p,e,q,N0,rdt,Nb,side,refraction,fshift);

        if(pex!=nullptr)
        {
        const auto t0 = kclk::now();
        pex->start(p,pgc,*e->N);
        // Phase 19: the halo exchange after a pipelined sweep (it waits for the ranks still in the sweep)
        if(pipelined_sweep())
        kst.exch += std::chrono::duration<double>(kclk::now()-t0).count();
        }
    }
}

void seastate_implicit::sweep(lexer *p, fdm_seastate *e, int q, const seastate_store *N0, double rdt,
                              const vector<float> &Nb, const int side[4], bool refraction, bool fshift)
{
    if(m1[q]<m0[q])
    return;

    // Phase 21b (A 767 1): the source terms once per iteration, evaluated by the sweep of the first quadrant (qfill) and
    // kept for the others; with the wind sea (wind input or DIA, the expensive terms), not in residual mode (FAS: the
    // residual of the present spectrum), not with the surfbeat model or the composite sweep of the mesh refinement
    // (its cells are solved one by one through solve(), whose terms are evaluated as before)
    suse = sonce && src!=nullptr && (src->param().wind || src->param().dia) && !second && rout==nullptr
           && vis==nullptr && visT==nullptr;
    if(suse)
    {
    scache_alloc(p,e);
    if(q==scache->qfill)
    ++sgen;
    }

    // Phase 14: pipelined across the ranks
    if(pipelined_sweep())
    {
    sweep_threads(p,e,q,N0,rdt,Nb,side,refraction,fshift);
    return;
    }

    // tiles (Phase 13, also with one thread: the upwind neighbours of a cell stay in the cache) and threads (A 798),
    // on level 0 and on the patches of the mesh refinement (covered cells skipped); in the composite sweep (A 797 1)
    // the visit of a covered cell gets the thread (visit_threads)
    // (not when the residuals of neighbouring cells are added into the same coarse cell, residual_add with s > 0)
    const int nbox = ranged ? (ri1-ri0+1)*(rj1-rj0+1) : p->knox*p->knoy;
    if((vis==nullptr || visT!=nullptr) && !second && !(rout!=nullptr && raccum && rshift>0)
       && nbox>=thread_min && (nthreads>1 || nbox>=4096))     // one thread: tiles from 16 tiles of 16 x 16 cells on
    {
    sweep_threads(p,e,q,N0,rdt,Nb,side,refraction,fshift);
    return;
    }

    const bool idown = (q==1 || q==2);
    const bool jdown = (q==2 || q==3);

    // the whole grid, or (mesh refinement) the interior of a patch
    const int ia = ranged ? ri0 : 0, ib = ranged ? ri1 : p->knox-1;
    const int ja = ranged ? rj0 : 0, jb = ranged ? rj1 : p->knoy-1;

    for(int ii=ia; ii<=ib; ++ii)
    for(int jj=ja; jj<=jb; ++jj)
    {
    i = idown ? ia+ib-ii : ii;
    j = jdown ? ja+jb-jj : jj;

        if(skp!=nullptr && (*skp)(i,j)==1)
        {
            // composite sweep (mesh refinement): the finer cells under this one, in its place
            if(vis!=nullptr)
            {
            const int si = i, sj = j;
            (*vis)(si,sj);
            i = si;
            j = sj;
            }
        continue;
        }

        if(e->wet(i,j)==1)
        {
            if(second)
            cell_surfbeat(p,e,q,N0,rdt,Nb,side,refraction);
            else
            cellw(p,e,q,i,j,N0!=nullptr ? N0->spec(i,j) : nullptr,rdt,Nb,side,refraction,fshift);
        }
    }
}

// threads (A 798): the cells with ii + jj = d (ii, jj in the order of the quadrant) depend only on
// the cells of the diagonals d-1 and d-2 (upwind neighbours, second-order fluxes) and are solved at
// the same time, each by a copy of the solver with its own source terms; the result is that of the
// serial sweep. The active ranges (spectral sparsity) are shared and set for all cells before, as
// the first visit of a cell sets them (and zeroes its bins below the threshold)
void seastate_implicit::sweep_threads(lexer *p, fdm_seastate *e, int q, const seastate_store *N0, double rdt,
                                      const vector<float> &Nb, const int side[4], bool refraction, bool fshift)
{
    const bool idown = (q==1 || q==2);
    const bool jdown = (q==2 || q==3);

    // the whole grid, or (mesh refinement) the interior of a patch
    const int ia = ranged ? ri0 : 0, ja = ranged ? rj0 : 0;
    const int ni = ranged ? ri1-ri0+1 : p->knox, nj = ranged ? rj1-rj0+1 : p->knoy;

    if(ni<=0 || nj<=0)
    return;

    const int nt = std::min(nthreads,std::max(ni,nj));
    const bool kba = pipelined_sweep();

    if(eps>0.0)
    {
    ranges_all(p,e,ia,ni,ja,nj,nt);

        // Phase 14: the first visits zero the bins below the threshold in all cells before the first sweep; the halo
        // takes them as the serial sweep sees them (once, on all ranks at the same sweep)
        if(kba && kfirst)
        kex->start(p,kpgc,*e->N);
    }
    kfirst = false;

    // copies of the solver and of the source terms for the threads 1..nt-1, made once (Phase 12) and given the
    // settings of this solver before every sweep
    worker_cache &W = workers;
    if(W.from!=src)
    {
    W.wk.clear();
    W.ws.clear();
    W.from = src;
    }
    while(int(W.wk.size())<nt)
    {
    W.wk.push_back(std::make_shared<seastate_implicit>(*this));
    W.ws.push_back(src!=nullptr ? std::make_shared<seastate_source>(*src) : nullptr);
    }

    std::vector<seastate_implicit*> sol(nt,this);

    for(int t=1; t<nt; ++t)
    {
    seastate_implicit *w = W.wk[t].get();
    w->copy_config(*this);
        if(src!=nullptr)
        {
        W.ws[t]->sync_param(*src);
        w->src = W.ws[t].get();
        }
    sol[t] = w;
    }

    // Phase 13: tiles of tb x tb cells (in the order of the quadrant) are the units of the threads. A tile waits for
    // its upwind tiles (the tile before it in i and the one before it in j); inside a tile the cells are solved in the
    // order of the serial sweep. So every cell has its four neighbours solved before or after it as in the serial
    // sweep, and the result is that of the serial sweep. Unlike the cell diagonals (Phase 7b), the tiles keep the
    // upwind neighbours of a cell in the cache of its thread: on the sparse swell spectra the cells took 22 % more time
    // with 2 threads than alone. The thread that finishes a tile goes on with the next tile in j if it is ready
    // Phase 19: pipelined, tiles of tbI x tbJ cells: thin across a rank boundary (pipeline)
    const int tbs = kba ? 16 : tile_size(ni,nj,nt);
    const int tbI = kba ? ktbI : tbs, tbJ = kba ? ktbJ : tbs;
    const int nI = (ni+tbI-1)/tbI, nJ = (nj+tbJ-1)/tbJ, ntile = nI*nJ;

    // Phase 14: the ranks upstream and downstream of this one in the order of the quadrant (-1: none). A tile at the
    // upstream boundary also waits for the halo cells of its rows (columns) from the rank upstream; a tile at the
    // downstream boundary sends its cells to the rank downstream when it is solved
    int upI = -1, dnI = -1, upJ = -1, dnJ = -1;
    if(kba)
    {
    upI = idown ? p->nb4 : p->nb1;
    dnI = idown ? p->nb1 : p->nb4;
    upJ = jdown ? p->nb2 : p->nb3;
    dnJ = jdown ? p->nb3 : p->nb2;
    upI = upI>=0 ? upI : -1; dnI = dnI>=0 ? dnI : -1;
    upJ = upJ>=0 ? upJ : -1; dnJ = dnJ>=0 ? dnJ : -1;
    }

    std::vector<std::atomic<int>> npred(ntile);
    for(int I=0; I<nI; ++I)
    for(int J=0; J<nJ; ++J)
    npred[size_t(I)*nJ+J].store((I>0 ? 1 : 0) + (J>0 ? 1 : 0) + (I==0 && upI>=0 ? 1 : 0) + (J==0 && upJ>=0 ? 1 : 0),
                                std::memory_order_relaxed);

    std::mutex qm;
    std::vector<int> ready;
    ready.reserve(ntile);
    if(npred[0].load(std::memory_order_relaxed)==0)
    ready.push_back(0);
    std::atomic<int> done(0);

    // Phase 14: messages of the boundary tiles. Axis 0: the rows of tile column J at the i boundary (cells ii = ni-1-n
    // downstream, the halo ii = -1-n upstream, n < L, in the order of the quadrant), axis 1: the columns of tile row I at
    // the j boundary. Cells without storage send zeros and ignore what they receive, as the halo exchange
    const int L = kba ? kex->halo_layers() : 0;
    const int nbin = e->grid->nbin;
    struct kmsg {int axis, idx; std::vector<float> buf; MPI_Request rq = MPI_REQUEST_NULL;};
    std::vector<kmsg> kin, kout;
    std::vector<int> outbox;            // kout entries packed by a thread, to be sent by thread 0
    size_t nposted = 0;

    auto cell_of = [&](int axis, int idx, int n, int k, bool halo, int &ci, int &cj)
    {
        const int ii = (axis==0) ? (halo ? -1-n : ni-1-n) : idx*tbI+k;
        const int jj = (axis==0) ? idx*tbJ+k : (halo ? -1-n : nj-1-n);
        ci = ia + (idown ? ni-1-ii : ii);
        cj = ja + (jdown ? nj-1-jj : jj);
    };
    auto length = [&](int axis, int idx)
    {
        return (axis==0) ? std::min(nj,(idx+1)*tbJ) - idx*tbJ : std::min(ni,(idx+1)*tbI) - idx*tbI;
    };
    // Phase 17: the sides of the rank downstream and upstream of it along an axis (0 i-, 1 i+, 2 j-, 3 j+)
    // (only in the composite sweep, which visits the covered cells; the FAS coarse sweeps of level 0 carry no patches)
    const int sdn[2] = {idown ? 0 : 1, jdown ? 2 : 3}, sup[2] = {idown ? 1 : 0, jdown ? 3 : 2};
    seastate_kba_extra *kx = (kba && visT!=nullptr) ? kext : nullptr;
    auto extra = [&](int axis, int idx, bool send) -> size_t
    {
        return kx!=nullptr ? kx->kba_size(q,send ? sdn[axis] : sup[axis],idx,axis==0 ? tbJ : tbI,send) : 0;
    };
    auto pack = [&](kmsg &m)
    {
        const int len = length(m.axis,m.idx);
        size_t c = 0;
        for(int n=0; n<L; ++n)
        for(int k=0; k<len; ++k)
        {
        int ci, cj;
        cell_of(m.axis,m.idx,n,k,false,ci,cj);
        const float *f = e->N->spec(ci,cj);
            if(f!=nullptr && e->N->active(ci,cj))
            std::copy(f,f+nbin,&m.buf[c]);
            else
            std::fill(&m.buf[c],&m.buf[c]+nbin,0.0f);
        c += nbin;
        }
        if(kx!=nullptr && c<m.buf.size())
        kx->kba_pack(q,sdn[m.axis],m.idx,m.axis==0 ? tbJ : tbI,&m.buf[c]);
    };
    auto unpack = [&](kmsg &m)
    {
        const int len = length(m.axis,m.idx);
        size_t c = 0;
        for(int n=0; n<L; ++n)
        for(int k=0; k<len; ++k)
        {
        int ci, cj;
        cell_of(m.axis,m.idx,n,k,true,ci,cj);
        float *f = e->N->spec(ci,cj);
        if(f!=nullptr && e->N->active(ci,cj))
        std::copy(&m.buf[c],&m.buf[c]+nbin,f);
        c += nbin;
        }
        if(kx!=nullptr && c<m.buf.size())
        kx->kba_unpack(q,sup[m.axis],m.idx,m.axis==0 ? tbJ : tbI,&m.buf[c]);
    };
    auto tag = [](int axis, int idx) {return 5000 + 2*idx + axis;};

    if(kba)
    {
        for(int axis=0; axis<2; ++axis)
        {
        const int up = axis==0 ? upI : upJ, dn = axis==0 ? dnI : dnJ, n = axis==0 ? nJ : nI;
            for(int idx=0; idx<n; ++idx)
            {
                if(up>=0)
                {
                kin.push_back(kmsg{axis,idx,std::vector<float>(size_t(L)*length(axis,idx)*nbin + extra(axis,idx,false))});
                MPI_Irecv(kin.back().buf.data(),int(kin.back().buf.size()),MPI_FLOAT,up,tag(axis,idx),kpgc->mpi_comm,&kin.back().rq);
                }
            if(dn>=0)
            kout.push_back(kmsg{axis,idx,std::vector<float>(size_t(L)*length(axis,idx)*nbin + extra(axis,idx,true))});
            }
        }
    outbox.reserve(kout.size());
    }
    // the kout entry of a boundary tile (axis, idx), -1: none
    const int nkoutI = (dnI>=0) ? nJ : 0;
    auto kout_of = [&](int axis, int idx) {return axis==0 ? (dnI>=0 ? idx : -1) : (dnJ>=0 ? nkoutI+idx : -1);};

    // Phase 14: the order of the tiles: towards the rank downstream in i the tile columns first (I fastest), towards
    // the one in j the tile rows (J fastest), towards both the diagonals I + J
    auto prio = [&](int k)
    {
        const long I = k/nJ, J = k%nJ;
        if(dnI>=0 && dnJ<0) return J*nI + I;
        if(dnJ>=0 && dnI<0) return I*nJ + J;
        return (I+J)*long(nI+nJ) + I;
    };

    // thread 0: send what the threads have packed, take in what has arrived (the tile it frees becomes ready)
    std::vector<int> kidx(kin.size());
    std::vector<MPI_Request> kreq(kin.size());
    auto pump = [&]()
    {
        std::vector<int> send;
        {
        std::lock_guard<std::mutex> lk(qm);
        send.swap(outbox);
        }
        for(int o : send)
        {
        kmsg &m = kout[o];
        MPI_Isend(m.buf.data(),int(m.buf.size()),MPI_FLOAT,m.axis==0 ? dnI : dnJ,tag(m.axis,m.idx),kpgc->mpi_comm,&m.rq);
        ++nposted;
        }

        // the messages in flight progress only within MPI calls (large messages wait for the receiver): test them
        for(kmsg &m : kout)
            if(m.rq!=MPI_REQUEST_NULL)
            {
            int f = 0;
            MPI_Test(&m.rq,&f,MPI_STATUS_IGNORE);
            }

        if(kin.empty())
        return;
        for(size_t r=0; r<kin.size(); ++r)
        kreq[r] = kin[r].rq;
        int nd = 0;
        MPI_Testsome(int(kin.size()),kreq.data(),&nd,kidx.data(),MPI_STATUSES_IGNORE);
        if(nd==MPI_UNDEFINED)
        return;
        for(int a=0; a<nd; ++a)
        {
        kmsg &m = kin[kidx[a]];
        m.rq = MPI_REQUEST_NULL;
        unpack(m);
        const int k = (m.axis==0) ? m.idx : m.idx*nJ;
            if(npred[k].fetch_sub(1,std::memory_order_acq_rel)==1)
            {
            std::lock_guard<std::mutex> lk(qm);
            ready.push_back(k);
            }
        }
    };

    // Phase 19: the time of the pipelined sweep per thread: the work of the tiles, the waiting for a ready tile before the
    // first one (the fill of the pipeline) and after it (the polling of the messages in it), the messages that thread 0
    // posts and tests after each of its tiles
    const auto ks0 = kclk::now();
    std::vector<double> tbusy(nt,0.0), tfill(nt,0.0), tidle(nt,0.0);
    double tcomm = 0.0;
    auto dsec = [](kclk::duration d) {return std::chrono::duration<double>(d).count();};

    auto work = [&](int t)
    {
        seastate_implicit *s = sol[t];
        int k = -1;
        auto mark = kclk::now();
        bool first = true;
        // the time since mark: waiting for a tile
        auto waited = [&]()
        {
            const auto n = kclk::now();
            (first ? tfill[t] : tidle[t]) += dsec(n-mark);
            mark = n;
        };
        auto pumpt = [&]()
        {
            waited();
            pump();
            const auto n = kclk::now();
            tcomm += dsec(n-mark);
            mark = n;
        };

        for(int spin=0;;)
        {
            if(k<0)
            {
                // while waiting for a tile: the polling of the messages counts as waiting
                if(kba && t==0)
                pump();

                {
                std::lock_guard<std::mutex> lk(qm);
                    if(!ready.empty())
                    {
                        if(kba)
                        {
                        // Phase 14: the ready tile that the ranks downstream need first
                        size_t b = 0;
                        for(size_t r=1; r<ready.size(); ++r)
                        if(prio(ready[r])<prio(ready[b]))
                        b = r;
                        std::swap(ready[b],ready.back());
                        }
                    k = ready.back();
                    ready.pop_back();
                    }
                }

                if(k<0)
                {
                    if(done.load(std::memory_order_acquire)>=ntile)
                    {
                    if(kba && t==0)
                    pump();
                    waited();
                    break;
                    }
                if(++spin>2000)
                std::this_thread::yield();
                continue;
                }
            spin = 0;
            waited();
            first = false;
            }

        const int I = k/nJ, J = k%nJ;
        const int ie = std::min(ni,(I+1)*tbI), je = std::min(nj,(J+1)*tbJ);
        const auto tw0 = kclk::now();

            for(int ii=I*tbI; ii<ie; ++ii)
            for(int jj=J*tbJ; jj<je; ++jj)
            {
            const int ci = ia + (idown ? ni-1-ii : ii);
            const int cj = ja + (jdown ? nj-1-jj : jj);

                if(skp!=nullptr && (*skp)(ci,cj)==1)
                {
                if(visT!=nullptr)
                (*visT)(t,ci,cj);
                }
                else if(e->wet(ci,cj)==1)
                s->cellw(p,e,q,ci,cj,N0!=nullptr ? N0->spec(ci,cj) : nullptr,rdt,Nb,side,refraction,fshift);
            }

        // Phase 14: the cells at the downstream boundary of the rank for the rank downstream
        if(kba)
        {
        const int oi = (I==nI-1) ? kout_of(0,J) : -1, oj = (J==nJ-1) ? kout_of(1,I) : -1;
        if(oi>=0) pack(kout[oi]);
        if(oj>=0) pack(kout[oj]);
            if(oi>=0 || oj>=0)
            {
            std::lock_guard<std::mutex> lk(qm);
            if(oi>=0) outbox.push_back(oi);
            if(oj>=0) outbox.push_back(oj);
            }
        }
        mark = kclk::now();
        tbusy[t] += dsec(mark-tw0);

        int next = -1;
        if(J+1<nJ && npred[k+1].fetch_sub(1,std::memory_order_acq_rel)==1)
        next = k+1;
            if(I+1<nI && npred[k+nJ].fetch_sub(1,std::memory_order_acq_rel)==1)
            {
                if(next<0)
                next = k+nJ;
                else
                {
                int other = k+nJ;
                if(kba && prio(other)<prio(next))
                std::swap(other,next);
                std::lock_guard<std::mutex> lk(qm);
                ready.push_back(other);
                }
            }
        done.fetch_add(1,std::memory_order_acq_rel);
        k = next;

        // Phase 14: thread 0 posts the messages after every tile
        if(kba && t==0)
        pumpt();
        }
    };

    const std::function<void(int)> job = work;
    seastate_pool::get().run(nt,job);

    if(kba)
    {
    const auto kf0 = kclk::now();
    if(nposted!=kout.size())
    pump();
        for(kmsg &m : kout)
        if(m.rq!=MPI_REQUEST_NULL)
        MPI_Wait(&m.rq,MPI_STATUS_IGNORE);
    const auto kf1 = kclk::now();

    ++kst.sweeps;
    kst.wall += dsec(kf1-ks0);
    kst.finish += dsec(kf1-kf0);
    kst.comm += tcomm;
        for(int t=0; t<nt; ++t)
        {
        kst.busy += tbusy[t];
        kst.fill += tfill[t];
        kst.idle += tidle[t];
        }
    }
}

bool seastate_implicit::pipelined_sweep() const
{
    if(kpgc==nullptr || kex==nullptr || ranged || second || rout!=nullptr)
    return false;

    // Phase 17: the composite sweep of the mesh refinement (covered cells visited with the index of the thread)
    if(skp!=nullptr || vis!=nullptr || visT!=nullptr)
    return kext!=nullptr && skp!=nullptr && visT!=nullptr;

    return true;
}

// Phase 14: the pipeline needs the same rows on both sides of an i boundary and the same columns on both sides of a j
// boundary (a Cartesian decomposition), and one tile size on all ranks
bool seastate_implicit::pipeline(lexer *p, ghostcell *pgc)
{
    int ok = 1;
    const int nb[4] = {p->nb1,p->nb4,p->nb3,p->nb2};
    for(int s=0; s<4; ++s)
    if(nb[s]>=0)
    {
    int mine[2] = {s<2 ? p->origin_j : p->origin_i, s<2 ? p->knoy : p->knox}, theirs[2] = {0,0};
    MPI_Sendrecv(mine,2,MPI_INT,nb[s],6100+s,theirs,2,MPI_INT,nb[s],6100+(s^1),pgc->mpi_comm,MPI_STATUS_IGNORE);
    if(theirs[0]!=mine[0] || theirs[1]!=mine[1])
    ok = 0;
    }
    int all = 0;
    MPI_Allreduce(&ok,&all,1,MPI_INT,MPI_MIN,pgc->mpi_comm);

    int tb = tile_size(p->knox,p->knoy,nthreads), tmin = tb;
    MPI_Allreduce(&tb,&tmin,1,MPI_INT,MPI_MIN,pgc->mpi_comm);
    tmin = std::max(tmin,1);

    // Phase 19: the tiles thin across the rank boundaries. A rank downstream in i starts when the rank upstream has
    // solved its first strip of tiles (the tile row J = 0 over all its columns), so with nJ tile rows the last of P ranks
    // in i starts (P - 1)/nJ of a sweep late, and the first one waits as long at the end of the sweep: nJ >= 8 (P - 1)
    // tile rows keep this below an eighth of a sweep. Along the sweep the tiles keep 16 cells (the cache of the upwind
    // neighbours). Ranks in i and j: the ranks with origin_j 0 (origin_i 0) of the Cartesian decomposition
    int cnt[2] = {p->origin_j==0 ? 1 : 0, p->origin_i==0 ? 1 : 0}, P[2] = {1,1};
    MPI_Allreduce(cnt,P,2,MPI_INT,MPI_SUM,pgc->mpi_comm);
    int nloc[2] = {p->knox,p->knoy}, nmin[2] = {1,1};
    MPI_Allreduce(nloc,nmin,2,MPI_INT,MPI_MIN,pgc->mpi_comm);
    int tI = tmin, tJ = tmin;
    if(P[0]>1)
    tJ = std::max(1,std::min(tmin,nmin[1]/(8*(P[0]-1))));
    if(P[1]>1)
    tI = std::max(1,std::min(tmin,nmin[0]/(8*(P[1]-1))));

    // validation 68: the tile sizes of the pipeline given as REEF3D_SEASTATE_TILES="tI,tJ" (the same on all ranks)
    if(const char *ev = std::getenv("REEF3D_SEASTATE_TILES"))
    {
    int a = 0, b = 0;
    if(std::sscanf(ev,"%d,%d",&a,&b)==2 && a>0 && b>0)
    {tI = a; tJ = b;}
    }

    if(all==1)
    {
    kpgc = pgc;
    ktbI = tI;
    ktbJ = tJ;
    }
    return all==1;
}

// tiles of the threaded sweep: 16 x 16 cells; with several threads smaller while the tile graph holds less than 3
// tiles per thread on average along its longest chain (nI nJ / (nI + nJ - 1))
int seastate_implicit::tile_size(int ni, int nj, int nt)
{
    int tb = 16;
        while(nt>1 && tb>2)
        {
        const double nI = std::ceil(double(ni)/tb), nJ = std::ceil(double(nj)/tb);
        if(nI*nJ/(nI+nJ-1.0)>=3.0*nt)
        break;
        tb /= 2;
        }
    return tb;
}

void seastate_implicit::run_threads(int nt, const std::function<void(int)> &f)
{
    seastate_pool::get().run(nt,f);
}

// the settings of another solver of the same grid (the threads of a sweep, Phase 12): everything the cell kernel
// reads besides its own work arrays
void seastate_implicit::copy_config(const seastate_implicit &o)
{
    for(int k=0; k<4; ++k)
    {
    Nbs[k] = o.Nbs[k];
    Nbs0[k] = o.Nbs0[k];
    }
    wU = o.wU;
    wD = o.wD;
    vN = o.vN;
    second = o.second;
    ranged = o.ranged;
    ri0 = o.ri0; ri1 = o.ri1; rj0 = o.rj0; rj1 = o.rj1;
    rwid = o.rwid;
    skp = o.skp;
    fcs = o.fcs;
    eps = o.eps;
    order2 = o.order2;
    rout = o.rout;
    raccum = o.raccum;
    rshift = o.rshift;
    ftau = o.ftau;
    ftil = o.ftil;
    fhold = o.fhold;
    nsrcit = o.nsrcit;
    srctol = o.srctol;
    vis = o.vis;
    wfld = o.wfld;
    visT = o.visT;
    dca = o.dca; dcax = o.dcax; dcay = o.dcay;
    pob = o.pob;
    nthreads = o.nthreads;
    thread_min = o.thread_min;
    sr = o.sr;
    slag = o.slag;
    lagon = o.lagon;
    sonce = o.sonce;
    suse = o.suse;
    sgen = o.sgen;
    scache = o.scache;
}

// Phase 21b (A 767 1): the slots of the wet cells and the layout of the kept rows, once
void seastate_implicit::scache_alloc(lexer *p, fdm_seastate *e)
{
    source_cache &c = *scache;
    if(c.ready)
    return;

    c.qfill = -1;
    c.wa = ndir;
    c.wb = -1;
    for(int q=0; q<4; ++q)
    if(m0[q]<=m1[q])
    {
    if(c.qfill<0)
    c.qfill = q;
    c.wa = std::min(c.wa,m0[q]);
    c.wb = std::max(c.wb,m1[q]);
    }

    c.width = 0;
    for(int q=0; q<4; ++q)
    {
    c.qoff[q] = -1;
        if(m0[q]<=m1[q] && q!=c.qfill)
        {
        c.qoff[q] = c.width;
        c.width += m1[q]-m0[q]+1;
        }
    }

    c.imin = p->imin; c.jmin = p->jmin; c.ni = p->imax; c.nj = p->jmax;
    c.slot.assign(size_t(c.ni)*c.nj,-1);

    int ns = 0;
    if(c.width>0)
    for(int ci=0; ci<p->knox; ++ci)
    for(int cj=0; cj<p->knoy; ++cj)
    if(e->wet(ci,cj)==1 && e->N->spec(ci,cj)!=nullptr)
    {
    const size_t x = size_t(ci-c.imin)*c.nj + (cj-c.jmin);
    if(x<c.slot.size())
    c.slot[x] = ns++;
    }

    c.PD.assign(size_t(ns)*2*nsig*c.width,0.0f);
    c.etot.assign(ns,0.0);
    c.tri.assign(ns,0);
    c.gen.assign(ns,-1);
    c.ready = true;
}

size_t seastate_implicit::source_cache_floats() const
{
    return scache->ready ? size_t(2*nsig*scache->width) : 0;
}

void seastate_implicit::preset_ranges(lexer *p, fdm_seastate *e)
{
    if(!(eps>0.0))
    return;

    const int ia = ranged ? ri0 : 0, ib = ranged ? ri1 : p->knox-1;
    const int ja = ranged ? rj0 : 0, jb = ranged ? rj1 : p->knoy-1;

    for(int ci=ia; ci<=ib; ++ci)
    for(int cj=ja; cj<=jb; ++cj)
    if(managed(p,e,ci,cj))
    ranges(p,e,ci,cj);
}

void seastate_implicit::solve(lexer *p, fdm_seastate *e, int q, int ci, int cj, const seastate_store *N0, double rdt,
                              const vector<float> &Nb, const int side[4], bool refraction, bool fshift)
{
    if(m1[q]<m0[q])
    return;

    i = ci;
    j = cj;

    // the composite sweep: the source terms evaluated in every sweep (A 767 is for the sweeps of a whole grid)
    suse = false;

    if(e->wet(i,j)==1)
    cellw(p,e,q,i,j,N0!=nullptr ? N0->spec(i,j) : nullptr,rdt,Nb,side,refraction,fshift);
}

// Phase 18: the solve of a cell, its time added into the work field
void seastate_implicit::cellw(lexer *p, fdm_seastate *e, int q, int ci, int cj, const float *N0c, double rdt,
                              const vector<float> &Nb, const int side[4], bool refraction, bool fshift)
{
    if(wfld==nullptr)
    {
    cell(p,e,q,ci,cj,N0c,rdt,Nb,side,refraction,fshift);
    return;
    }

    const auto t0 = std::chrono::steady_clock::now();
    cell(p,e,q,ci,cj,N0c,rdt,Nb,side,refraction,fshift);
    (*wfld)(ci,cj) += std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
}

void seastate_implicit::ranges_get(lexer *p, int ci, int cj, double *v)
{
    const int n = ranges_size();
    if(n==0)
    return;
    if(rg.empty() || ci<rimin || ci>=rimin+rni || cj<rjmin || cj>=rjmin+rnj)
    {
    std::fill(v,v+n,0.0);
    return;
    }
    const size_t c = size_t(ci-rimin)*rnj + (cj-rjmin);
    v[0] = double(rv[c]);
    for(int k=0; k<4; ++k)
    v[1+k] = double(rb[c*4+k]);
    for(int k=0; k<4*nsig; ++k)
    v[5+k] = double(rg[c*4*nsig+k]);
}

void seastate_implicit::ranges_put(lexer *p, int ci, int cj, const double *v)
{
    if(ranges_size()==0)
    return;
    ranges_alloc(p);
    if(ci<rimin || ci>=rimin+rni || cj<rjmin || cj>=rjmin+rnj)
    return;
    const size_t c = size_t(ci-rimin)*rnj + (cj-rjmin);
    rv[c] = uint16_t(v[0]);
    for(int k=0; k<4; ++k)
    rb[c*4+k] = uint16_t(v[1+k]);
    for(int k=0; k<4*nsig; ++k)
    rg[c*4*nsig+k] = uint16_t(v[5+k]);
}

// Phase 23 (A 768 1): ecWAM's wave stress and u* of the last time step per cell
double seastate_implicit::lag_key = 0.0;

void seastate_implicit::lagging(bool on, lexer *p)
{
    lagon = on;
    shared_lag &L = *slag;
    if(!on || !L.tw.empty())
    return;
    L.ni = p->imax; L.nj = p->jmax; L.imin = p->imin; L.jmin = p->jmin;
    const size_t n = size_t(L.ni)*L.nj;
    L.tw.assign(n,0.0); L.twd.assign(n,0.0); L.us.assign(n,0.0);
    L.nw.assign(n,0.0); L.nwd.assign(n,0.0); L.nus.assign(n,0.0);
    L.key.assign(n,-1.0e300);
}

long seastate_implicit::lag_cell(lexer *p, int ci, int cj)
{
    shared_lag &L = *slag;
    if(L.tw.empty())
    return -1;
    if(ci<L.imin || ci>=L.imin+L.ni || cj<L.jmin || cj>=L.jmin+L.nj)
    return -1;
    const long c = long(ci-L.imin)*L.nj + (cj-L.jmin);
    // the first evaluation of the cell in a new time step: the state the last step left
    if(L.key[c]!=lag_key)
    {
        if(L.key[c]>-1.0e300)
        {
        L.tw[c] = L.nw[c];
        L.twd[c] = L.nwd[c];
        L.us[c] = L.nus[c];
        }
    L.key[c] = lag_key;
    }
    return c;
}

void seastate_implicit::lag_start(lexer *p, int ci, int cj)
{
    if(!lagon || src==nullptr || !src->lags())
    return;
    const long c = lag_cell(p,ci,cj);
    if(c>=0)
    src->lag(slag->tw[c],slag->twd[c],slag->us[c]);
}

void seastate_implicit::lag_end(lexer *p, int ci, int cj)
{
    if(!lagon || src==nullptr || !src->lags())
    return;
    const long c = lag_cell(p,ci,cj);
    if(c<0)
    return;
    shared_lag &L = *slag;
    L.nw[c] = src->lag_tauw();
    L.nwd[c] = src->lag_tauwdir();
    L.nus[c] = src->lag_ustar();
}

// checkpoints: the state the last completed step left (the step in progress has not begun at a checkpoint)
void seastate_implicit::lag_get(lexer *p, int ci, int cj, double *v)
{
    if(lag_size()==0)
    return;
    v[0] = v[1] = v[2] = 0.0;
    const shared_lag &L = *slag;
    if(L.tw.empty() || ci<L.imin || ci>=L.imin+L.ni || cj<L.jmin || cj>=L.jmin+L.nj)
    return;
    const size_t c = size_t(ci-L.imin)*L.nj + (cj-L.jmin);
    if(L.key[c]<=-1.0e300)
    return;
    v[0] = L.nw[c];
    v[1] = L.nwd[c];
    v[2] = L.nus[c];
}

void seastate_implicit::lag_put(lexer *p, int ci, int cj, const double *v)
{
    if(lag_size()==0)
    return;
    shared_lag &L = *slag;
    if(L.tw.empty() || ci<L.imin || ci>=L.imin+L.ni || cj<L.jmin || cj>=L.jmin+L.nj)
    return;
    const size_t c = size_t(ci-L.imin)*L.nj + (cj-L.jmin);
    // as written by the last step: the next evaluation (a new step) takes it over
    L.nw[c] = v[0];
    L.nwd[c] = v[1];
    L.nus[c] = v[2];
    L.key[c] = -1.0e299;
}

namespace
{
    using neighbour = seastate_neighbour;

    // fs: side of the neighbour that faces the cell (mesh refinement: spectrum on the face of a covered cell)
    void set_neighbour(lexer *p, fdm_seastate *e, int ni, int nj, const vector<float> &Nb, const vector<float> *const Nbs[4], const int side[4], neighbour &nb,
                       sliceint *skp = nullptr, const seastate_faces *fcs = nullptr, int fs = 0)
    {
        nb = neighbour();

        if(e->wet(ni,nj)==1)
        {
        nb.N  = e->N->spec(ni,nj);

            if(fcs!=nullptr && skp!=nullptr && (*skp)(ni,nj)==1)
            {
            const float *f = fcs->face(ni,nj,fs);
            if(f!=nullptr)
            nb.N = f;
            }

        nb.cg = e->cg->spec(ni,nj);
        nb.U  = e->U(ni,nj);
        nb.V  = e->V(ni,nj);
        return;
        }

        const int gi = ni + p->origin_i;
        const int gj = nj + p->origin_j;

        int s=-1;
        if(gi<0)
        s=0;
        else if(gi>=p->gknox)
        s=1;
        else if(gj<0)
        s=2;
        else if(gj>=p->gknoy)
        s=3;

        if(s>=0 && side[s]==1 && !Nb.empty())
        nb.N = Nb.data();

        // surfbeat (x-) and boundary series (A 711 3): spectrum of the boundary cell
        if(s>=0 && side[s]==1 && Nbs[s]!=nullptr)
        {
        const int k = (s<2) ? nj : ni;
            if(k>=0 && k<((s<2) ? p->knoy : p->knox))
            nb.N = Nbs[s]->data() + size_t(k)*size_t(e->N->nbins());
        }

        if(s>=0 && side[s]==2)
        nb.self = true;
    }

    inline double pos(double a) {return a>0.0 ? a : 0.0;}
    inline double neg(double a) {return a<0.0 ? a : 0.0;}

    // second-order correction of the upwind flux through a face with velocity c:
    // upwind cell u, the cell upstream of it uu, the downstream cell dn (active cells, else 0)
    inline double tvd(double c, const float *uu, const float *u, const float *dn, int b)
    {
        if(uu==nullptr || u==nullptr || dn==nullptr)
        return 0.0;

        const double du = double(dn[b]) - double(u[b]);

        if(!(std::fabs(du)>1.0e-30))
        return 0.0;

        const double r = (double(u[b]) - double(uu[b]))/du;
        const double phi = (r + std::fabs(r))/(1.0 + std::fabs(r));

        return 0.5*c*phi*du;
    }

    inline const float *active(fdm_seastate *e, int i, int j)
    {
        return e->wet(i,j)==1 ? e->N->spec(i,j) : nullptr;
    }
}

// ------------------------------------------------------------------------------------------------
// spectral sparsity (A 795): per cell, frequency and quadrant the range of directions that hold
// energy (rg, lo<<8 | hi, 0xFFFF empty); bins outside the range are zero. A sweep solves per
// frequency the window of the bins that can receive energy: the own range (one more direction on
// each side with refraction, the next frequencies with frequency shift, the lower frequencies of the
// triads), the ranges of the four neighbours, the ends of the quadrant next to energy in the
// neighbouring quadrants; the whole quadrant where wind input or DIA can fill any bin.

// total energy of a spectrum
double seastate_implicit::energy(const seastate_grid &g, const float *N) const
{
    double et = 0.0;
    for(int l=0; l<nsig; ++l)
    {
    double s = 0.0;
    const float *Nl = N + g.bin(l,0);
    for(int m=0; m<ndir; ++m)
    s += double(Nl[m])*g.wth[m];
    et += s*g.sig[l]*g.dsig[l]*g.dtheta;
    }
    return et;
}

// the cell holds its active ranges: solved by this solver (interior or patch box), wet, not covered
bool seastate_implicit::managed(lexer *p, fdm_seastate *e, int ci, int cj) const
{
    const int ia = ranged ? ri0 : 0, ib = ranged ? ri1 : p->knox-1;
    const int ja = ranged ? rj0 : 0, jb = ranged ? rj1 : p->knoy-1;

    if(ci<ia || ci>ib || cj<ja || cj>jb)
    return false;

    if(e->wet(ci,cj)!=1 || e->N->spec(ci,cj)==nullptr)
    return false;

    return skp==nullptr || (*skp)(ci,cj)==0;
}

void seastate_implicit::ranges_alloc(lexer *p)
{
    if(!rg.empty())
    return;

    rni = p->imax; rnj = p->jmax; rimin = p->imin; rjmin = p->jmin;
    rg.assign(size_t(rni)*rnj*4*nsig,0xFFFF);
    rv.assign(size_t(rni)*rnj,uint16_t(0));
    rb.assign(size_t(rni)*rnj*4,uint16_t(0xFFFF));
}

// the first visits of the cells of a box (ranges of the present spectrum), shared out by rows to nt threads: each
// writes only the ranges and the spectrum of its own cells
void seastate_implicit::ranges_all(lexer *p, fdm_seastate *e, int ia, int ni, int ja, int nj, int nt)
{
    ranges_alloc(p);

    const std::function<void(int)> job = [&](int t)
    {
        for(int ci=ia+t; ci<ia+ni; ci+=nt)
        for(int cj=ja; cj<ja+nj; ++cj)
        if(rv[size_t(ci-rimin)*rnj + (cj-rjmin)]==0 && managed(p,e,ci,cj))
        ranges(p,e,ci,cj);
    };
    seastate_pool::get().run(std::max(1,std::min(nt,ni)),job);
}

uint16_t *seastate_implicit::ranges(lexer *p, fdm_seastate *e, int ci, int cj)
{
    ranges_alloc(p);

    const size_t c = size_t(ci-rimin)*rnj + (cj-rjmin);
    uint16_t *r = &rg[c*4*nsig];

    if(rv[c]==0)
    {
    // first visit: ranges of the present spectrum (all quadrants), zero below the threshold
    const seastate_grid &g = *e->grid;
    float *N = e->N->spec(ci,cj);
    const double et = energy(g,N);

        for(int qq=0; qq<4; ++qq)
        {
        const int ma = m0[qq], nq = m1[qq]-m0[qq]+1;
            for(int l=0; l<nsig; ++l)
            {
            float *Nl = N + g.bin(l,ma);
            const double t = eps*et/(g.sig[l]*g.dsig[l]*g.dtheta);
            int lo = nq, hi = -1;
                for(int n=0; n<nq; ++n)
                if(double(Nl[n])>t)
                {
                lo = std::min(lo,n);
                hi = n;
                }
                for(int n=0; n<nq; ++n)
                if(n<lo || n>hi)
                Nl[n] = 0.0f;
            r[qq*nsig+l] = (hi<0) ? 0xFFFF : uint16_t((lo<<8) | hi);
            }
        }

    rv[c] = 1 | 2;
    for(int qq=0; qq<4; ++qq)
    summary(ci,cj,qq,r);
    }

    return r;
}

// quadrant summary of a cell: energy in the quadrant, in its first and its last direction
void seastate_implicit::summary(int ci, int cj, int q, const uint16_t *r)
{
    const size_t c = size_t(ci-rimin)*rnj + (cj-rjmin);
    const int nq = m1[q]-m0[q]+1;
    bool any=false, lo=false, hi=false;
    int la=nsig, lb=-1;

    for(int l=0; l<nsig; ++l)
    {
    const uint16_t rr = r[q*nsig+l];
    if(rr==0xFFFF)
    continue;
    any = true;
    la = std::min(la,l);
    lb = l;
    if((rr>>8)==0) lo = true;
    if(int(rr&0xFF)==nq-1) hi = true;
    }

    uint16_t v = rv[c] & uint16_t(~((4<<q) | (64<<(2*q)) | (128<<(2*q))));
    if(any) v |= uint16_t(4<<q);
    if(lo)  v |= uint16_t(64<<(2*q));
    if(hi)  v |= uint16_t(128<<(2*q));
    rv[c] = v;
    rb[c*4+q] = any ? uint16_t((la<<8) | lb) : uint16_t(0xFFFF);
}

static inline void hull(int &lo, int &hi, uint16_t r)
{
    if(r==0xFFFF)
    return;
    lo = std::min(lo,int(r>>8));
    hi = std::max(hi,int(r&0xFF));
}

void seastate_implicit::windows(lexer *p, fdm_seastate *e, int q, int ci, int cj, const float *N,
                                const neighbour &W, const neighbour &E, const neighbour &S, const neighbour &Nn,
                                bool rf, bool fs, bool all)
{
    const seastate_grid &g = *e->grid;
    const int ma = m0[q], nq = m1[q]-m0[q]+1;
    const uint16_t *r = ranges(p,e,ci,cj);

    // the whole quadrant where the source terms can fill any bin, or reflection (obstacles, coasts) can
    // turn energy of any quadrant into this one
    if(all || (src!=nullptr && src->fills_spectrum()))
    {
    for(int l=0; l<nsig; ++l)
    {
    wlo[l] = 0;
    whi[l] = nq-1;
    }
    return;
    }

    const int qm = (q+3)%4, qp = (q+1)%4;
    const int nqm = m1[qm]-m0[qm]+1;
    const uint16_t sc = rv[size_t(ci-rimin)*rnj + (cj-rjmin)];
    const int tr = (src!=nullptr && (sc & 2)) ? src->triad_reach() : 0;

    // neighbours: their ranges if this solver holds them, else the positive bins of their spectrum
    const neighbour *nb[4] = {&W,&E,&S,&Nn};
    const int nbi[4] = {ci-1,ci+1,ci,ci}, nbj[4] = {cj,cj,cj-1,cj+1};
    const uint16_t *nr[4] = {nullptr,nullptr,nullptr,nullptr};

    bool maybe = (sc & (4<<q)) || (rf && ((sc & (128<<(2*qm))) || (sc & (64<<(2*qp)))));

    for(int k=0; k<4; ++k)
    if(nb[k]->N!=nullptr)
    {
        if(nb[k]->cg!=nullptr && managed(p,e,nbi[k],nbj[k]) && nb[k]->N==e->N->spec(nbi[k],nbj[k]))
        {
        nr[k] = ranges(p,e,nbi[k],nbj[k]);
        if(rv[size_t(nbi[k]-rimin)*rnj + (nbj[k]-rjmin)] & (4<<q))
        maybe = true;
        }
        else
        maybe = true;
    }

    for(int l=0; l<nsig; ++l)
    {
    wlo[l] = 0;
    whi[l] = -1;
    }

    // nothing in this quadrant, its neighbours or (refraction) at the ends of the next quadrants
    if(!maybe)
    return;

    // band of frequencies that can receive energy
    const size_t cc = size_t(ci-rimin)*rnj + (cj-rjmin);
    int La = nsig, Lb = -1;
    auto band = [&](uint16_t b, int dn, int up)
    {
        if(b==0xFFFF)
        return;
        La = std::min(La,std::max(int(b>>8)-dn,0));
        Lb = std::max(Lb,std::min(int(b&0xFF)+up,nsig-1));
    };

    band(rb[cc*4+q],fs ? 1 : 0,(fs ? 1 : 0) - tr);
    if(rf)
    {
    if(sc & (128<<(2*qm))) band(rb[cc*4+qm],0,0);
    if(sc & (64<<(2*qp)))  band(rb[cc*4+qp],0,0);
    }
    for(int k=0; k<4; ++k)
    if(nb[k]->N!=nullptr)
    {
        if(nr[k]!=nullptr)
        band(rb[(size_t(nbi[k]-rimin)*rnj + (nbj[k]-rjmin))*4+q],0,0);
        else
        {
        La = 0;
        Lb = nsig-1;
        }
    }

    for(int l=La; l<=Lb; ++l)
    {
    int lo = nq, hi = -1;

    const uint16_t own = r[q*nsig+l];
    hull(lo,hi,own);

        // refraction: one more direction on each side, and the quadrant ends next to energy in the
        // neighbouring quadrants
        if(rf)
        {
            if(own!=0xFFFF)
            {
            lo = std::max(lo-1,0);
            hi = std::min(hi+1,nq-1);
            }

        const uint16_t rm = r[qm*nsig+l], rp = r[qp*nsig+l];
        if(rm!=0xFFFF && int(rm&0xFF)==nqm-1)
        {lo = std::min(lo,0); hi = std::max(hi,0);}
        if(rp!=0xFFFF && int(rp>>8)==0)
        {lo = std::min(lo,nq-1); hi = std::max(hi,nq-1);}
        }

        // frequency shift: the next frequencies
        if(fs)
        {
        if(l>0)      hull(lo,hi,r[q*nsig+l-1]);
        if(l<nsig-1) hull(lo,hi,r[q*nsig+l+1]);
        }

        // triads: the lower frequencies that feed this one
        for(int k=l+tr; k<l; ++k)
        if(k>=0)
        hull(lo,hi,r[q*nsig+k]);

        // inflow from the neighbours
        for(int k=0; k<4; ++k)
        if(nb[k]->N!=nullptr)
        {
            if(nr[k]!=nullptr)
            hull(lo,hi,nr[k][q*nsig+l]);
            else
            {
            const float *Nl = nb[k]->N + g.bin(l,ma);
            int n=0;
            while(n<nq && !(Nl[n]>0.0f)) ++n;
            if(n<nq)
            {
            int m=nq-1;
            while(!(Nl[m]>0.0f)) --m;
            lo = std::min(lo,n);
            hi = std::max(hi,m);
            }
            }
        }

    wlo[l] = lo;
    whi[l] = hi;
    }
}

// after the solve: the new ranges of the solved windows, zero below the threshold at their ends
void seastate_implicit::keep(lexer *p, fdm_seastate *e, int q, int ci, int cj)
{
    const seastate_grid &g = *e->grid;
    const int ma = m0[q];
    uint16_t *r = ranges(p,e,ci,cj);
    float *N = e->N->spec(ci,cj);

    for(int l=lmin; l<=lmax; ++l)
    {
    const int na = wlo[l], nb = whi[l];
    if(na>nb)
    {
    r[q*nsig+l] = 0xFFFF;
    continue;
    }

    float *Nl = N + g.bin(l,ma);
    int lo = nb+1, hi = na-1;

        for(int n=na; n<=nb; ++n)
        if(double(Nl[n])>thr[l])
        {
        if(lo>nb) lo = n;
        hi = n;
        }

        for(int n=na; n<=nb; ++n)
        if(n<lo || n>hi)
        Nl[n] = 0.0f;

    r[q*nsig+l] = (hi<lo) ? 0xFFFF : uint16_t((lo<<8) | hi);
    }

    summary(ci,cj,q,r);
}

// second-order fluxes: spectrum of an active cell within the solved area or the ring around it,
// not covered by a finer grid
const float *seastate_implicit::usable(lexer *p, fdm_seastate *e, int ci, int cj) const
{
    const int ia = ranged ? ri0-rwid : -2, ib = ranged ? ri1+rwid : p->knox+1;
    const int ja = ranged ? rj0-rwid : -2, jb = ranged ? rj1+rwid : p->knoy+1;

    if(ci<ia || ci>ib || cj<ja || cj>jb)
    return nullptr;

    if(e->wet(ci,cj)!=1)
    return nullptr;

    if(skp!=nullptr && (*skp)(ci,cj)==1)
    return nullptr;

    return e->N->spec(ci,cj);
}

// solution of frequency l (directions wlo..whi of the quadrant starting at ma) into N: zero after
// a vanishing pivot, the action density limiter of the source terms (|change per iteration| <= dmax,
// dmax < 0: none), clipped at zero with the second-order fluxes; with spectral sparsity the bins below
// the threshold at the ends of the window are set to zero (the active range of the frequency)
void seastate_implicit::store(float *N, int l, int ma, int nq, double dmax, bool clip)
{
    const int na = wlo[l], nb = whi[l], o = l*ndir;
    float *Nl = N + l*ndir + ma;

    if(na>nb)
    return;

    if(tok[l]==0)
    {
    for(int n=na; n<=nb; ++n)
    Nl[n] = 0.0f;
    return;
    }

    // the limiter bounds the change of the iteration: with source iterations (A 738) the change from the
    // spectrum of the cell before its first solve, not the change of each solve
    const float *Nb = (Nlim!=nullptr) ? Nlim + l*ndir + ma : Nl;

    if(dmax>=0.0)
    for(int n=na; n<=nb; ++n)
    {
    const double v = double(Nb[n]);
    tso[o+n] = std::min(std::max(tso[o+n],v-dmax),v+dmax);
    }

    if(clip)
    for(int n=na; n<=nb; ++n)
    tso[o+n] = std::max(tso[o+n],0.0);

    if(nup!=nullptr)
    for(int n=na; n<=nb; ++n)
    tso[o+n] = std::min(tso[o+n],nup[l]);

    for(int n=na; n<=nb; ++n)
    Nl[n] = float(tso[o+n]);
}

void seastate_implicit::cell(lexer *p, fdm_seastate *e, int q, int ci, int cj, const float *N0, double rdt,
                             const vector<float> &Nb, const int side[4], bool refraction, bool fshift)
{
    // the cell (ci,cj); local i, j (IP, JP), the solver works in several threads (A 798)
    int i = ci, j = cj;

    const seastate_grid &g = *e->grid;

    float *N = e->N->spec(i,j);
    const float *kc = e->kw->spec(i,j);
    const float *cgc = e->cg->spec(i,j);

    const int ic = i, jc = j;

    neighbour W,E,S,Nn;
    set_neighbour(p,e,ic-1,jc,Nb,Nbs,side,W,skp,fcs,1);
    set_neighbour(p,e,ic+1,jc,Nb,Nbs,side,E,skp,fcs,0);
    set_neighbour(p,e,ic,jc-1,Nb,Nbs,side,S,skp,fcs,3);
    set_neighbour(p,e,ic,jc+1,Nb,Nbs,side,Nn,skp,fcs,2);

    const double rdx = 1.0/p->DXN[IP];
    const double rdy = 1.0/p->DYN[JP];
    const double rdth = 1.0/g.dtheta;

    // obstacles (A 722) on the faces of the cell: the inflow through a face is transmitted with Kt^2
    const seastate_obstacle::face *oW = nullptr, *oE = nullptr, *oS = nullptr, *oN = nullptr;
    if(pob!=nullptr)
    {
    oW = pob->east(ic-1,jc);
    oE = pob->east(ic,jc);
    oS = pob->north(ic,jc-1);
    oN = pob->north(ic,jc);
    if(oW) {W.tf = oW->kt2; W.tff = pob->kt2f(oW);}
    if(oE) {E.tf = oE->kt2; E.tff = pob->kt2f(oE);}
    if(oS) {S.tf = oS->kt2; S.tff = pob->kt2f(oS);}
    if(oN) {Nn.tf = oN->kt2; Nn.tff = pob->kt2f(oN);}
    }
    const double rdxW = rdx*W.tf, rdxE = rdx*E.tf, rdyS = rdy*S.tf, rdyN = rdy*Nn.tf;
    auto reflecting = [&](const seastate_obstacle::face *f) {return f!=nullptr && (f->kr2>0.0f || f->fq>=0);};
    const bool refl = reflecting(oW) || reflecting(oE) || reflecting(oS) || reflecting(oN);

    // diffraction (A 718): Ca and its gradient of the cell and Ca of the neighbours
    const float *cac = (dca!=nullptr) ? dca->spec(ic,jc) : nullptr;
    const float *cax = (dcax!=nullptr) ? dcax->spec(ic,jc) : nullptr;
    const float *cay = (dcay!=nullptr) ? dcay->spec(ic,jc) : nullptr;
    const bool dif = (cac!=nullptr && cax!=nullptr && cay!=nullptr);
    if(dif)
    {
    if(W.cg)  W.ca  = dca->spec(ic-1,jc);
    if(E.cg)  E.ca  = dca->spec(ic+1,jc);
    if(S.cg)  S.ca  = dca->spec(ic,jc-1);
    if(Nn.cg) Nn.ca = dca->spec(ic,jc+1);
    }

    const double d = e->depth(ic,jc);
    const double U = e->U(ic,jc), V = e->V(ic,jc);
    const bool kin = (e->refr(ic,jc)==1);

    const double ddx = e->ddx(ic,jc), ddy = e->ddy(ic,jc);
    const double dUdx = e->dUdx(ic,jc), dUdy = e->dUdy(ic,jc), dVdx = e->dVdx(ic,jc), dVdy = e->dVdy(ic,jc);
    const double dddt = e->dddt(ic,jc);

    const int ma = m0[q], mb = m1[q], nq = mb-ma+1;

    // frequency shift only where c_sigma can be nonzero (depth change in time or along the current,
    // current gradients); without it the frequencies are independent
    const double dep = dddt + U*ddx + V*ddy;
    const bool fs = fshift && kin && (dep!=0.0 || dUdx!=0.0 || dUdy!=0.0 || dVdx!=0.0 || dVdy!=0.0);
    const bool rf = refraction && kin;

    // wind field: the wind of the cell
    if(src!=nullptr && wU!=nullptr)
    src->set_wind((*wU)(ic,jc),(*wD)(ic,jc));

    // vegetation field: the stems of the cell
    if(src!=nullptr && vN!=nullptr)
    src->set_vegetation((*vN)(ic,jc));

    // directions solved per frequency: the whole quadrant, or (spectral sparsity, A 795) the window
    // of the bins that can hold energy
    const bool sparse = (eps>0.0);
    int wmin = 0, wmax = nq-1;

    if(sparse)
    {
    windows(p,e,q,ic,jc,N,W,E,S,Nn,rf || dif,fs,refl);

    wmin = nq; wmax = -1;
    lmin = nsig; lmax = -1;
        for(int l=0; l<nsig; ++l)
        if(wlo[l]<=whi[l])
        {
        wmin = std::min(wmin,wlo[l]);
        wmax = std::max(wmax,whi[l]);
        lmin = std::min(lmin,l);
        lmax = l;
        }

        // nothing can arrive in this quadrant
        if(wmax<0)
        {
        i = ic;
        j = jc;
        return;
        }
    }
    else
    {
    lmin = 0;
    lmax = nsig-1;
        for(int l=0; l<nsig; ++l)
        {
        wlo[l] = 0;
        whi[l] = nq-1;
        }
    }

    // source iterations per cell (A 738): the source terms are evaluated again from the spectrum just
    // solved and the cell is solved again
    // (residual mode, FAS: one evaluation)
    const int nsi = (rout!=nullptr) ? 1 : nsrcit;
    double esit = (nsi>1) ? energy(g,N) : 0.0, dsit = 0.0;

    if(nsi>1)
    Nsi0.assign(N,N+g.nbin);
    Nlim = (nsi>1) ? Nsi0.data() : nullptr;

    for(int sit=0; sit<nsi; ++sit)
    {
    if(nsi>1 && sit==1)
    Nsi1.assign(N,N+g.nbin);

    // Phase 21b (A 767 1): the first evaluation in a cell either evaluates the terms for the iteration (the sweep of
    // qfill: sfill) or takes them (the other quadrants: skept, if the sweep of qfill evaluated them in this iteration)
    int cslot = -1;
    bool sfill = false, skept = false;
    double setot = 0.0;
    bool stri = false;

    if(src!=nullptr && suse && sit==0)
    {
    cslot = scache->slot_of(ic,jc);
        if(cslot>=0)
        {
        sfill = (q==scache->qfill);
        skept = !sfill && scache->gen[cslot]==sgen;
        }
    }

    // maximum energy (A 737 1, SWAN SINTGRL), then the source terms from the latest spectrum,
    // for the directions that are solved; with spectral sparsity the row sums over the active ranges
    if(src!=nullptr)
    {
        if(sparse)
        {
        const uint16_t *r = ranges(p,e,ic,jc);
            for(int l=0; l<nsig; ++l)
            {
            double s = 0.0;
                for(int qq=0; qq<4; ++qq)
                {
                const uint16_t rr = r[qq*nsig+l];
                if(rr==0xFFFF)
                continue;
                const float *Nl = N + g.bin(l,m0[qq]);
                const double *wq = g.wth.data() + m0[qq];
                for(int n=int(rr>>8); n<=int(rr&0xFF); ++n)
                s += double(Nl[n])*wq[n];
                }
            rsum[l] = s;
            }
        src->set_rows(rsum.data());
        }

    const bool capped = src->cap(N,d);

        // Phase 21b (A 767 1): the terms kept from the sweep of qfill in this iteration, for the whole quadrant (the
        // window may still grow, widen); a cell scaled down to the maximum energy just now is evaluated again
        if(skept && !capped)
        {
        const source_cache &c = *scache;
        const size_t w = size_t(c.width);
        const float *cP = &c.PD[size_t(cslot)*2*nsig*w], *cD = cP + size_t(nsig)*w;
        const int o = c.qoff[q];

            for(int l=0; l<nsig; ++l)
            {
            const float *pl = cP + size_t(l)*w + o, *dl = cD + size_t(l)*w + o;
            double *Pl = P.data() + g.bin(l,ma), *Dl = D.data() + g.bin(l,ma);
                for(int n=0; n<nq; ++n)
                {
                Pl[n] = double(pl[n]);
                Dl[n] = double(dl[n]);
                }
            }

        src->limits(kc,cgc);
        src->drop_rows();
        setot = c.etot[cslot];
        stri = (c.tri[cslot]!=0);
        }
        else
        {
        skept = false;

            // the sweep of qfill: all directions and frequencies, the bins of the other quadrants kept
            if(sfill)
            {
            source_cache &c = *scache;
            lag_start(p,ic,jc);
            src->compute(N,d,kc,cgc,P.data(),D.data(),c.wa,c.wb,0,nsig-1);
            lag_end(p,ic,jc);

            const size_t w = size_t(c.width);
            float *cP = &c.PD[size_t(cslot)*2*nsig*w], *cD = cP + size_t(nsig)*w;

                for(int qq=0; qq<4; ++qq)
                if(c.qoff[qq]>=0)
                {
                const int mq = m0[qq], nqq = m1[qq]-m0[qq]+1, o = c.qoff[qq];
                    for(int l=0; l<nsig; ++l)
                    {
                    const double *Pl = P.data() + g.bin(l,mq), *Dl = D.data() + g.bin(l,mq);
                    float *pl = cP + size_t(l)*w + o, *dl = cD + size_t(l)*w + o;
                        for(int n=0; n<nqq; ++n)
                        {
                        pl[n] = float(Pl[n]);
                        dl[n] = float(Dl[n]);
                        }
                    }
                }

            c.etot[cslot] = src->Etot;
            c.tri[cslot] = src->triads_active() ? 1 : 0;
            c.gen[cslot] = sgen;
            }
            else
            {
            lag_start(p,ic,jc);
            src->compute(N,d,kc,cgc,P.data(),D.data(),ma+wmin,ma+wmax,lmin,lmax);
            lag_end(p,ic,jc);
            }

        setot = src->Etot;
        stri = src->triads_active();
        }
    }

    // FAS coarse-grid correction (A 758): the source -tau of the coarse cell, its loss part proportional to
    // N (exact at the restricted spectrum, where the coarse problem starts); Phase 15, with the spectral sparsity
    // (A 795) of the finer grid: a bin that is zero in all children is held at zero. The finer grid keeps the bins
    // below its threshold empty; filled on the coarse grid (wind input, DIA) they changed the coarse spectra at the
    // converged fine solution, the corrections did not vanish and grew after about 100 iterations
    if(ftau!=nullptr)
    {
    const float *T = ftau->spec(ic,jc), *Nt = (ftil!=nullptr) ? ftil->spec(ic,jc) : nullptr;

        if(src==nullptr)
        for(int l=lmin; l<=lmax; ++l)
        for(int n=wlo[l]; n<=whi[l]; ++n)
        P[g.bin(l,ma+n)] = D[g.bin(l,ma+n)] = 0.0;

        if(T!=nullptr && Nt!=nullptr)
        for(int l=lmin; l<=lmax; ++l)
        for(int n=wlo[l]; n<=whi[l]; ++n)
        {
        const int b = g.bin(l,ma+n);
        const double s = -double(T[b]);
            if(fhold && !(Nt[b]>0.0f))
            D[b] += 1.0e20;
            else if(s>0.0)
            P[b] += s;
            else if(s<0.0 && Nt[b]>1.0e-30f)
            D[b] -= s/double(Nt[b]);
        }
    }

    const double *lim = (src!=nullptr) ? src->limit() : nullptr;

    // ST4 (A 732 3): the upper bound of the spectrum per frequency (ecWAM FLMAX, Phillips alpha 0.031)
    nup = (src!=nullptr) ? src->upper() : nullptr;

    // triads of this cell (the window of the next sweeps includes the lower frequencies that feed a bin)
    if(sparse && src!=nullptr)
    {
    const size_t c = size_t(ic-rimin)*rnj + (jc-rjmin);
    rv[c] = stri ? (rv[c] | 2) : (rv[c] & ~2);
    }

    // threshold of the sparse spectrum: bins with E < A 795 E_tot of the cell are set to zero
    if(sparse)
    {
    const double et = (src!=nullptr) ? setot : energy(g,N);
    for(int l=0; l<nsig; ++l)
    thr[l] = eps*et/(g.sig[l]*g.dsig[l]*g.dtheta);
    }

    // coefficients of the cell that do not depend on the bin: per direction of the quadrant (current
    // shear, refraction factors of the theta faces) and per frequency (depth refraction
    // sig/sinh(2kd), c_sigma); then one pass over the bins

    for(int n=0; n<nq; ++n)
    {
    const int m = ma+n;
    const double cs = g.costh[m], sn = g.sinth[m];

    qc[n] = cs;
    qs[n] = sn;
    cur[n] = cs*cs*dUdx + sn*cs*(dUdy+dVdx) + sn*sn*dVdy;
    rth[n] = rdth/g.wth[m];

        if(!rf)
        tAp[n] = tCp[n] = tAm[n] = tCm[n] = 0.0;
        else
        {
        const double sp = g.sinthf[m],  cp_ = g.costhf[m];
        const double sm = g.sinthl[m],  cm  = g.costhl[m];

        tAp[n] = sp*ddx - cp_*ddy;
        tCp[n] = sp*cp_*(dUdx-dVdy) + sp*sp*dVdx - cp_*cp_*dUdy;
        tAm[n] = sm*ddx - cm*ddy;
        tCm[n] = sm*cm*(dUdx-dVdy) + sm*sm*dVdx - cm*cm*dUdy;
        }
    }

    if(kin)
    for(int l=0; l<nsig; ++l)
    if(fs || wlo[l]<=whi[l])
    Ar[l] = seastate_refraction(g.sig[l],kc[l],d);

    // diffraction: c_theta = Ca (depth refraction) + current refraction + cg dCa/dn, n = (-sin, cos) the
    // normal to the direction (as SWAN), at the faces of the directions of the quadrant
    if(dif)
    for(int l=0; l<nsig; ++l)
    if(wlo[l]<=whi[l])
    {
    const double A = (kin ? Ar[l] : 0.0)*double(cac[l]);
    const double ax = double(cgc[l])*double(cax[l]), ay = double(cgc[l])*double(cay[l]);
    double *dp = &cdp[size_t(l)*ndir], *dm = &cdm[size_t(l)*ndir];

        for(int n=0; n<nq; ++n)
        {
        const int m = ma+n;
        dp[n] = A*tAp[n] + tCp[n] - ax*g.sinthf[m]  + ay*g.costhf[m];
        dm[n] = A*tAm[n] + tCm[n] - ax*g.sinthl[m]  + ay*g.costhl[m];
        }
    }

    // c_sigma at the cell centre, per frequency and direction of the quadrant (zero without
    // frequency shift)
    if(!fs)
    for(int l=0; l<nsig; ++l)
    for(int n=0; n<nq; ++n)
    csg[size_t(l)*ndir+n] = 0.0;
    else
    for(int l=0; l<nsig; ++l)
    {
    const double a = double(kc[l])*Ar[l]*dep;
    const double bk = double(cgc[l])*double(kc[l]);
    double *c = &csg[size_t(l)*ndir];

        for(int n=0; n<nq; ++n)
        c[n] = a - bk*cur[n];
    }

    // neighbour spectra: zero spectrum where there is no inflow
    const float *NW = W.N ? W.N : zero.data();
    const float *NE = E.N ? E.N : zero.data();
    const float *NS = S.N ? S.N : zero.data();
    const float *NN = Nn.N ? Nn.N : zero.data();
    const bool self = W.self || E.self || S.self || Nn.self;
    const float *Nold = (N0!=nullptr) ? N0 : N;
    const double *Dsrc = (src!=nullptr || ftau!=nullptr) ? D.data() : zerod.data();
    const double *Psrc = (src!=nullptr || ftau!=nullptr) ? P.data() : zerod.data();

    // second-order geographic fluxes (A 796 2): the cells up to two away (active, solved or the ring
    // around a patch, not covered by a finer grid)
    const bool so = (order2 && second==false);
    const float *aW=nullptr, *aE=nullptr, *aS=nullptr, *aN=nullptr, *aWW=nullptr, *aEE=nullptr, *aSS=nullptr, *aNN=nullptr;

    if(so)
    {
    aW  = usable(p,e,ic-1,jc); aE  = usable(p,e,ic+1,jc); aS  = usable(p,e,ic,jc-1); aN  = usable(p,e,ic,jc+1);
    aWW = usable(p,e,ic-2,jc); aEE = usable(p,e,ic+2,jc); aSS = usable(p,e,ic,jc-2); aNN = usable(p,e,ic,jc+2);
    }

    // with obstacles and coasts (A 722, A 723) the faces next to a blocked face stay first order: the
    // upwind cells of the outflow faces (o) and of the inflow faces (i) without a blocked face between them
    const float *aWo = aW, *aEo = aE, *aSo = aS, *aNo = aN;
    const float *aWi = aW, *aEi = aE, *aSi = aS, *aNi = aN;
    if(so && pob!=nullptr)
    {
    const bool bW = oW!=nullptr, bE = oE!=nullptr, bS = oS!=nullptr, bN = oN!=nullptr;
    const bool bWW = pob->east(ic-2,jc)!=nullptr, bEE = pob->east(ic+1,jc)!=nullptr;
    const bool bSS = pob->north(ic,jc-2)!=nullptr, bNN = pob->north(ic,jc+1)!=nullptr;
    if(bW || bE) {aWo = nullptr; aEo = nullptr;}
    if(bS || bN) {aSo = nullptr; aNo = nullptr;}
    if(bW || bWW) aWi = nullptr;
    if(bE || bEE) aEi = nullptr;
    if(bS || bSS) aSi = nullptr;
    if(bN || bNN) aNi = nullptr;
    }

    // A: the tridiagonal systems in theta of all frequencies (arrays [l*ndir + n], directions
    // wlo..whi of each frequency); the inflow from the next lower frequency (c_sigma > 0) is added
    // in C with its latest value
    auto assemble = [&](int l, int na, int nb)
    {
    const double cgl = cgc[l];
    const double A = kin ? Ar[l] : 0.0;
    const double rdsig = 1.0/g.dsig[l];

    // geographic face velocities: (c_cell + c_neighbour)/2 = cg_face cos(theta) + U_face; with
    // diffraction Ca cg
    double cgW = W.cg  ? 0.5*(cgl + double(W.cg[l]))  : cgl;
    double cgE = E.cg  ? 0.5*(cgl + double(E.cg[l]))  : cgl;
    double cgS = S.cg  ? 0.5*(cgl + double(S.cg[l]))  : cgl;
    double cgN = Nn.cg ? 0.5*(cgl + double(Nn.cg[l])) : cgl;
    const double UW = W.cg  ? 0.5*(U + W.U)  : U, UE = E.cg  ? 0.5*(U + E.U)  : U;
    const double VS = S.cg  ? 0.5*(V + S.V)  : V, VN = Nn.cg ? 0.5*(V + Nn.V) : V;

    if(dif)
    {
    const double cdl = double(cac[l])*cgl;
    auto face = [&](const neighbour &nb) {return nb.cg ? 0.5*(cdl + (nb.ca ? double(nb.ca[l]) : 1.0)*double(nb.cg[l])) : cdl;};
    cgW = face(W);
    cgE = face(E);
    cgS = face(S);
    cgN = face(Nn);
    }

    const double *cc = &csg[size_t(l)*ndir];
    const double *cm = (l>0)      ? &csg[size_t(l-1)*ndir] : cc;
    const double *cu = (l<nsig-1) ? &csg[size_t(l+1)*ndir] : cc;

    const int b0 = g.bin(l,ma);
    const float *Nlp = (l<nsig-1) ? N + g.bin(l+1,ma) : zero.data();
    const int o = l*ndir;

    const double *__restrict Dl = Dsrc + b0, *__restrict Pl = Psrc + b0;
    const float *__restrict Nol = Nold + b0, *__restrict NWl = NW + b0, *__restrict NEl = NE + b0;
    const float *__restrict NSl = NS + b0, *__restrict NNl = NN + b0, *__restrict Nlp_ = Nlp;
    const double *__restrict qc_ = qc.data(), *__restrict qs_ = qs.data(), *__restrict rth_ = rth.data();
    const double *__restrict tAp_ = tAp.data(), *__restrict tCp_ = tCp.data(), *__restrict tAm_ = tAm.data(), *__restrict tCm_ = tCm.data();
    const double *__restrict cc_ = cc, *__restrict cu_ = cu, *__restrict cm_ = cm;
    const double *__restrict cdp_ = &cdp[size_t(l)*ndir], *__restrict cdm_ = &cdm[size_t(l)*ndir];
    double *__restrict di_ = &tdi[o], *__restrict la_ = &tla[o], *__restrict up_ = &tup[o];

    // inflow coefficients of the faces, with the transmission of structures per frequency (A 725)
    const double rxW = W.tff ? rdx*double(W.tff[l]) : rdxW, rxE = E.tff ? rdx*double(E.tff[l]) : rdxE;
    const double ryS = S.tff ? rdy*double(S.tff[l]) : rdyS, ryN = Nn.tff ? rdy*double(Nn.tff[l]) : rdyN;
    double *__restrict rh_ = &trh[o], *__restrict si_ = &tsi[o];

        // branch-free, vectorised over the directions
        #pragma GCC ivdep
        for(int n=na; n<=nb; ++n)
        {
        const double cs = qc_[n], sn = qs_[n];

        // c_sigma at the faces l+1/2 and l-1/2 (zero without frequency shift; range ends: outflow only)
        const double csc = cc_[n];
        const double csu = 0.5*(csc + cu_[n]);
        const double csl = 0.5*(csc + cm_[n]);

        // c_theta at the faces m+1/2 and m-1/2 (zero without refraction)
        const double ctp = dif ? cdp_[n] : A*tAp_[n] + tCp_[n];
        const double ctm = dif ? cdm_[n] : A*tAm_[n] + tCm_[n];

        const double cxw = cgW*cs + UW, cxe = cgE*cs + UE;
        const double cys = cgS*sn + VS, cyn = cgN*sn + VN;

        // diagonal: 1/dt + outflow (+ D)
        di_[n] = rdt + (pos(cxe) - neg(cxw))*rdx + (pos(cyn) - neg(cys))*rdy
                     + (pos(csu) - neg(csl))*rdsig + Dl[n] + (pos(ctp) - neg(ctm))*rth_[n];

        // right-hand side: old time level + inflow from the neighbours and from l+1 (+ P)
        rh_[n] = rdt*double(Nol[n]) + Pl[n]
               + pos(cxw)*rxW*double(NWl[n]) - neg(cxe)*rxE*double(NEl[n])
               + pos(cys)*ryS*double(NSl[n]) - neg(cyn)*ryN*double(NNl[n])
               - neg(csu)*rdsig*double(Nlp_[n]);

        // theta: within the window implicit (tridiagonal)
        la_[n] = -pos(ctm)*rth_[n];
        up_[n] =  neg(ctp)*rth_[n];
        si_[n] = pos(csl)*rdsig;
        }

    // theta: outside the window (other quadrants, empty bins) with the latest values
    {
    const double ctm0 = dif ? cdm_[na] : A*tAm[na] + tCm[na], ctp1 = dif ? cdp_[nb] : A*tAp[nb] + tCp[nb];
    const int mlo = g.prev(ma+na), mhi = g.next(ma+nb);       // -1: the end of a compact sector (no inflow)
    if(mlo>=0)
    rh_[na] += pos(ctm0)*rth[na]*double(N[g.bin(l,mlo)]);
    if(mhi>=0)
    rh_[nb] -= neg(ctp1)*rth[nb]*double(N[g.bin(l,mhi)]);
    la_[na] = 0.0;
    up_[nb] = 0.0;
    }

        // obstacles (A 722): the flux leaving the cell through a reflecting face in the direction
        // theta_i = 2 alpha - theta re-enters it in theta (specular at the obstacle line), Kr^2 of it;
        // N(theta_i) with the latest values, linear between the bin centres
        if(refl)
        {
        const double pi2 = 6.28318530717958647692;
        const float *Nl0 = N + g.bin(l,0);

        auto ninterp = [&](double t)
        {
            t = std::fmod(t,pi2);
            if(t<0.0)
            t += pi2;
            // compact sector (A 717): directions from theta0 on, nothing outside the sector
            if(!g.periodic)
            {
            double d = std::fmod(t-g.theta0,pi2);
            if(d<0.0)
            d += pi2;
            const double x = d/g.dtheta;
            if(x>double(ndir-1))
            return 0.0;
            const int a = std::min(int(x),ndir-2);
            const double w = x-double(a);
            return (1.0-w)*double(Nl0[a]) + w*double(Nl0[a+1]);
            }
            const int b = int(std::upper_bound(g.theta.begin(),g.theta.end(),t) - g.theta.begin());
            const int m1 = (b==0) ? ndir-1 : b-1, m2 = (b==ndir) ? 0 : b;
            double d1 = g.theta[m1], d2 = g.theta[m2];
            if(d1>t) d1 -= pi2;
            if(d2<t) d2 += pi2;
            const double w = (d2>d1) ? (t-d1)/(d2-d1) : 0.0;
            return (1.0-w)*double(Nl0[m1]) + w*double(Nl0[m2]);
        };

        // the reflection works on copies of the values it shares with the vectorised loop above: the
        // lambda is not inlined, and locals captured by reference would live in memory, which keeps
        // that loop from being vectorised
        const int ra = na, rb = nb, rl = l;
        double *const rh = rh_;
        const double fcW = cgW, fcE = cgE, fcS = cgS, fcN = cgN, fuW = UW, fuE = UE, fvS = VS, fvN = VN;

        auto reflect = [&,ra,rb,rl,rh,fcW,fcE,fcS,fcN,fuW,fuE,fvS,fvN](const seastate_obstacle::face *f, int side)
        {
            if(!reflecting(f))
            return;

            const double kr2 = pob->kr2(f,rl);
            if(!(kr2>0.0))
            return;

            // diffuse reflection (A 726, A 723 pown): the incident directions around the specular one
            int nd = 0;
            const float *wd = pob->diffuse(f,nd);

            auto outflow = [&](double ti)
            {
                if(side==0) return -(fcW*std::cos(ti) + fuW);
                if(side==1) return  (fcE*std::cos(ti) + fuE);
                if(side==2) return -(fcS*std::sin(ti) + fvS);
                return (fcN*std::sin(ti) + fvN);
            };

            const double r = (side<2) ? rdx : rdy;

            for(int n=ra; n<=rb; ++n)
            {
            const double tn = g.theta[ma+n], ti = 2.0*double(f->alpha) - tn;

                // the reflected direction leaves the line on the incident side; on the staircase of faces it
                // may also cross a face of the cell, where it is transmitted and reflected again
                if(wd==nullptr)
                {
                const double cout = outflow(ti);
                if(cout>0.0)
                rh[n] += kr2*cout*r*ninterp(ti);
                }
                else
                {
                double fl = 0.0;
                    for(int k=-nd; k<=nd; ++k)
                    {
                    const double tk = ti + k*g.dtheta, cout = outflow(tk);
                    if(cout>0.0)
                    fl += double(wd[nd+k])*cout*ninterp(tk);
                    }
                rh[n] += kr2*r*fl;
                }
            }
        };

        reflect(oW,0);
        reflect(oE,1);
        reflect(oS,2);
        reflect(oN,3);
        }

        // zero-gradient sides: inflow of the cell's own spectrum, implicit while the non-theta part
        // of the diagonal stays at least half of its value (M-matrix), otherwise with the latest value
        if(self)
        for(int n=na; n<=nb; ++n)
        {
        const double cs = qc[n], sn = qs[n];
        const double cxw = cgW*cs + UW, cxe = cgE*cs + UE;
        const double cys = cgS*sn + VS, cyn = cgN*sn + VN;
        const double ctp = dif ? cdp_[n] : A*tAp[n] + tCp[n], ctm = dif ? cdm_[n] : A*tAm[n] + tCm[n];
        const double dg = di_[n] - (pos(ctp) - neg(ctm))*rth[n];

        double aself = 0.0;
        if(W.self)  aself += pos(cxw)*rdx;
        if(E.self)  aself -= neg(cxe)*rdx;
        if(S.self)  aself += pos(cys)*rdy;
        if(Nn.self) aself -= neg(cyn)*rdy;

            if(aself>0.0)
            {
                if(dg-aself>=0.5*dg)
                di_[n] -= aself;
                else
                rh_[n] += aself*double(N[b0+n]);
            }
        }

        // second-order upwind geographic fluxes (A 796 2, as SWAN SORDUP): the face value is extrapolated
        // from the two upwind cells, F = c (3/2 N_up - 1/2 N_upup). Implicit in the cell, the upwind cells
        // are solved before it in the sweep; faces without two active upwind cells stay first order. Not
        // monotone: the solution is clipped at zero (store)
        if(so)
        for(int n=na; n<=nb; ++n)
        {
        const int b = b0+n;
        const double cs = qc[n], sn = qs[n];
        const double cxw = cgW*cs + UW, cxe = cgE*cs + UE;
        const double cys = cgS*sn + VS, cyn = cgN*sn + VN;
        double dd = 0.0, rr = 0.0;

        // x+ face
        if(cxe>0.0 && aWo)         {dd += 0.5*cxe*rdx;  rr += 0.5*cxe*rdx*double(aWo[b]);}
        if(cxe<0.0 && aEi && aEE)  {rr -= 0.5*cxe*rdx*(double(aEi[b]) - double(aEE[b]));}
        // x- face
        if(cxw>0.0 && aWi && aWW)  {rr += 0.5*cxw*rdx*(double(aWi[b]) - double(aWW[b]));}
        if(cxw<0.0 && aEo)         {dd -= 0.5*cxw*rdx;  rr -= 0.5*cxw*rdx*double(aEo[b]);}
        // y+ face
        if(cyn>0.0 && aSo)         {dd += 0.5*cyn*rdy;  rr += 0.5*cyn*rdy*double(aSo[b]);}
        if(cyn<0.0 && aNi && aNN)  {rr -= 0.5*cyn*rdy*(double(aNi[b]) - double(aNN[b]));}
        // y- face
        if(cys>0.0 && aSi && aSS)  {rr += 0.5*cys*rdy*(double(aSi[b]) - double(aSS[b]));}
        if(cys<0.0 && aNo)         {dd -= 0.5*cys*rdy;  rr -= 0.5*cys*rdy*double(aNo[b]);}

        di_[n] += dd;
        rh_[n] += rr;
        }
    };

    for(int l=lmin; l<=lmax; ++l)
    if(wlo[l]<=whi[l])
    assemble(l,wlo[l],whi[l]);

    // residual mode (FAS, A 758): r = b - A N of the cell equations at the latest spectra (the pseudo
    // time step drops out: N_old = N), for the directions of the quadrant; no solve
    // (residual_add, Phase 12: added into the cell (ic >> s, jc >> s) of R, no zeroing: the sum over the children of a
    // coarse cell of the uniform-grid FAS)
    if(rout!=nullptr)
    {
    float *R = raccum ? rout->spec(ic>>rshift,jc>>rshift) : rout->spec(ic,jc);
    const int ltail = (src!=nullptr) ? src->tail_from() : nsig;

        if(R!=nullptr)
        {
            if(!raccum)
            for(int l=0; l<nsig; ++l)
            for(int n=0; n<nq; ++n)
            R[g.bin(l,ma+n)] = 0.0f;

            for(int l=lmin; l<=lmax; ++l)
            {
            const int o = l*ndir, na = wlo[l], nb = whi[l];
            const float *Nl = N + g.bin(l,ma);
            const float *Nlm = (l>0) ? N + g.bin(l-1,ma) : nullptr;

                for(int n=na; n<=nb; ++n)
                {
                double r = trh[o+n] - tdi[o+n]*double(Nl[n]);
                if(Nlm!=nullptr) r += tsi[o+n]*double(Nlm[n]);
                if(n>na) r -= tla[o+n]*double(Nl[n-1]);
                if(n<nb) r -= tup[o+n]*double(Nl[n+1]);
                // Phase 15: with the spectral sparsity a zero bin is inactive (below the threshold, set to zero by
                // the sweep), its equation does not hold at the converged solution: no residual
                if(sparse && !(Nl[n]>0.0f))
                r = 0.0;
                // Phase 20: the bins of the diagnostic tail (A 732 3) are set after the solve, their equations do not hold either
                if(l>=ltail)
                r = 0.0;
                if(raccum)
                R[g.bin(l,ma+n)] += float(r);
                else
                R[g.bin(l,ma+n)] = float(r);
                }
            }
        }

    i = ic;
    j = jc;
    return;
    }

    // B: elimination factors of all frequencies (independent, the inner loop over l hides the
    // latency of the divisions); a frequency with a vanishing pivot is set to zero
    // (with spectral sparsity the windows are short: one frequency after the other, in C)
    if(!sparse)
    {
    for(int l=lmin; l<=lmax; ++l)
    tok[l] = 1;

    for(int n=wmin; n<=wmax; ++n)
    for(int l=lmin; l<=lmax; ++l)
    {
    if(n<wlo[l] || n>whi[l])
    continue;

    const int x = l*ndir + n;
    const double den = tdi[x] - (n>wlo[l] ? tla[x]*tcp[x-1] : 0.0);

        if(!(den>1.0e-300))
        {
        tok[l] = 0;
        trd[x] = tcp[x] = 0.0;
        continue;
        }

    trd[x] = 1.0/den;
    tcp[x] = tup[x]*trd[x];
    }
    }

    // one frequency on its own: factors, the inflow from l-1 (frequency shift), substitution
    auto solve1 = [&](int l)
    {
        const int o = l*ndir, na = wlo[l], nb = whi[l];
        if(na>nb)
        return;

        tok[l] = 1;
        for(int n=na; n<=nb; ++n)
        {
        const int x = o+n;
        const double den = tdi[x] - (n>na ? tla[x]*tcp[x-1] : 0.0);
            if(!(den>1.0e-300))
            {
            tok[l] = 0;
            return;
            }
        trd[x] = 1.0/den;
        tcp[x] = tup[x]*trd[x];
        }

        if(fs && l>0)
        {
        const float *Nlm = N + g.bin(l-1,ma);
        for(int n=na; n<=nb; ++n)
        trh[o+n] += tsi[o+n]*double(Nlm[n]);
        }

        for(int n=na; n<=nb; ++n)
        tdp[o+n] = (trh[o+n] - (n>na ? tla[o+n]*tdp[o+n-1] : 0.0))*trd[o+n];

        tso[o+nb] = tdp[o+nb];
        for(int n=nb-1; n>=na; --n)
        tso[o+n] = tdp[o+n] - tcp[o+n]*tso[o+n+1];
    };

    // spectral sparsity: where the solution at an end of the window is above the threshold (energy
    // refracted further than the window within this solve), the window grows and the frequency is
    // solved again, with the source terms of the new bins
    auto widen = [&](int l)
    {
        for(int tries=0; tries<8; ++tries)
        {
        const int na = wlo[l], nb = whi[l];
        if(na>nb || tok[l]==0)
        return;

        // estimate of the bin beyond each end: its inflow through the theta face from the end bin
        // over the diagonal of the end bin
        const float *Nl = N + g.bin(l,ma);
        const double A = kin ? Ar[l] : 0.0;
        const double ctm0 = dif ? cdm[size_t(l)*ndir+na] : A*tAm[na] + tCm[na];
        const double ctp1 = dif ? cdp[size_t(l)*ndir+nb] : A*tAp[nb] + tCp[nb];
        const int o = l*ndir;
        const bool lo = (na>0    && -neg(ctm0)*rth[na-1]*double(Nl[na])>thr[l]*tdi[o+na]);
        const bool hi = (nb<nq-1 &&  pos(ctp1)*rth[nb+1]*double(Nl[nb])>thr[l]*tdi[o+nb]);
        if(!lo && !hi)
        return;

        const int w = std::max(2,nb-na+1);
        const int na2 = lo ? std::max(0,na-w) : na, nb2 = hi ? std::min(nq-1,nb+w) : nb;

            // (A 767 1: the kept terms hold the whole quadrant)
            if(src!=nullptr && !skept)
            {
            // ST4 (A 732 3): u*, the growth rates and the cut-off of the compute of this cell
            src->set_rows(rsum.data());
            src->freeze(true);
            src->compute(N,d,kc,cgc,P.data(),D.data(),ma+na2,ma+nb2,l,l);
            src->freeze(false);
            }

        wlo[l] = na2;
        whi[l] = nb2;
        wmin = std::min(wmin,na2);
        wmax = std::max(wmax,nb2);
        assemble(l,na2,nb2);
        solve1(l);
        store(N,l,ma,nq,lim!=nullptr ? lim[l] : -1.0,so);
        }
    };

    // C: forward and back substitution. Without frequency shift the frequencies are independent;
    // with it, the inflow from l-1 takes the value just solved (Gauss-Seidel in sigma)
    if(sparse)
    for(int l=lmin; l<=lmax; ++l)
    {
    solve1(l);
    store(N,l,ma,nq,lim!=nullptr ? lim[l] : -1.0,so);
    widen(l);
    }
    else if(!fs)
    {
        for(int n=wmin; n<=wmax; ++n)
        for(int l=lmin; l<=lmax; ++l)
        {
        if(n<wlo[l] || n>whi[l])
        continue;

        const int x = l*ndir + n;
        tdp[x] = (trh[x] - (n>wlo[l] ? tla[x]*tdp[x-1] : 0.0))*trd[x];
        }

        for(int n=wmax; n>=wmin; --n)
        for(int l=lmin; l<=lmax; ++l)
        {
        if(n<wlo[l] || n>whi[l])
        continue;

        const int x = l*ndir + n;
        tso[x] = (n<whi[l]) ? tdp[x] - tcp[x]*tso[x+1] : tdp[x];
        }

        for(int l=lmin; l<=lmax; ++l)
        store(N,l,ma,nq,lim!=nullptr ? lim[l] : -1.0,so);
    }
    else
    for(int l=lmin; l<=lmax; ++l)
    {
    const int o = l*ndir, na = wlo[l], nb = whi[l];

        if(na<=nb)
        {
            if(l>0)
            {
            const float *Nlm = N + g.bin(l-1,ma);
            for(int n=na; n<=nb; ++n)
            trh[o+n] += tsi[o+n]*double(Nlm[n]);
            }

            for(int n=na; n<=nb; ++n)
            tdp[o+n] = (trh[o+n] - (n>na ? tla[o+n]*tdp[o+n-1] : 0.0))*trd[o+n];

            tso[o+nb] = tdp[o+nb];
            for(int n=nb-1; n>=na; --n)
            tso[o+n] = tdp[o+n] - tcp[o+n]*tso[o+n+1];
        }

    store(N,l,ma,nq,lim!=nullptr ? lim[l] : -1.0,so);
    }

        // ST4 (A 732 3): the diagnostic tail above the cut-off of the source terms of this solve
        if(src!=nullptr)
        src->project(N,kc,cgc,ma,mb);

        // the cell has reached the balance of its source terms with the inflow: the distance of its
        // energy to the fixed point, estimated from the change d of the last solve and the contraction
        // rho = d/d_previous as d/(1-rho), relative below A 738
        if(nsi>1)
        {
        const double en = energy(g,N), d = std::fabs(en-esit);

        // no contraction (the sources and the solve of the cell oscillate or stall, e.g. steep breaking
        // dissipation, sources against the limiter): the spectrum of the first solve, end of the source
        // iterations (a fixed point of the outer iteration is then a solution of the cell, as without
        // source iterations; the mean of the last two spectra used before held two-cycles of the cell)
        if(sit>0 && d>=0.95*dsit && d>srctol*std::max(en,1.0e-30))
        {
        std::copy(Nsi1.begin(),Nsi1.end(),N);
        break;
        }

        const double rho = (sit>0 && dsit>0.0) ? std::min(d/dsit,0.95) : 0.0;
        if(d<=(1.0-rho)*srctol*std::max(en,1.0e-30))
        break;
        esit = en;
        dsit = d;
        }
    }

    Nlim = nullptr;

    // new active ranges of the cell
    if(sparse)
    keep(p,e,q,ic,jc);

    i = ic;
    j = jc;
}

/*--------------------------------------------------------------------
Surfbeat (A 775 2): one frequency, Crank-Nicolson in time for the
transport (theta = 1/2: no numerical diffusion from the time stepping,
the wave groups keep their shape over long distances), second-order
geographic fluxes (van Leer, deferred correction), the source terms
implicit (Patankar). Old values: N0 (the spectra at the start of the
step, halo included) and the boundary rows of the previous step.
Positive for c dt/dx <= 2 (the host's time step is limited by the
long-wave celerity, which exceeds c_g); N is clipped at zero.
--------------------------------------------------------------------*/

void seastate_implicit::cell_surfbeat(lexer *p, fdm_seastate *e, int q, const seastate_store *N0s, double rdt,
                                      const vector<float> &Nb, const int side[4], bool refraction)
{
    const seastate_grid &g = *e->grid;
    const double th = 0.5, th1 = 1.0-th;
    const int nbin = g.nbin;

    float *N = e->N->spec(i,j);
    const float *kc = e->kw->spec(i,j);
    const float *cgc = e->cg->spec(i,j);

    const int ic = i, jc = j;
    const float *Nc0 = (N0s!=nullptr) ? N0s->spec(ic,jc) : N;

    neighbour W,E,S,Nn;
    set_neighbour(p,e,ic-1,jc,Nb,Nbs,side,W);
    set_neighbour(p,e,ic+1,jc,Nb,Nbs,side,E);
    set_neighbour(p,e,ic,jc-1,Nb,Nbs,side,S);
    set_neighbour(p,e,ic,jc+1,Nb,Nbs,side,Nn);

    // old values of the neighbours: active cells from N0, the x- rows of the previous step, Nb
    auto old = [&](const neighbour &nb, int ni, int nj) -> const float*
    {
        if(nb.N==nullptr)
        return nullptr;

        if(nb.cg!=nullptr)
        return (N0s!=nullptr) ? N0s->spec(ni,nj) : nb.N;

        if(ni+p->origin_i<0 && Nbs0[0]!=nullptr && nj>=0 && nj<p->knoy)
        return Nbs0[0]->data() + size_t(nj)*size_t(nbin);

        return nb.N;
    };

    const float *W0 = old(W,ic-1,jc), *E0 = old(E,ic+1,jc), *S0 = old(S,ic,jc-1), *N0n = old(Nn,ic,jc+1);

    // second order: active neighbours and active cells two cells away, new and old values
    const float *aW = W.cg ? W.N : nullptr, *aE = E.cg ? E.N : nullptr, *aS = S.cg ? S.N : nullptr, *aN = Nn.cg ? Nn.N : nullptr;
    const float *oW = W.cg ? W0 : nullptr, *oE = E.cg ? E0 : nullptr, *oS = S.cg ? S0 : nullptr, *oN = Nn.cg ? N0n : nullptr;
    const float *aWW = active(e,ic-2,jc), *aEE = active(e,ic+2,jc), *aSS = active(e,ic,jc-2), *aNN = active(e,ic,jc+2);
    const float *oWW = (aWW && N0s) ? N0s->spec(ic-2,jc) : aWW, *oEE = (aEE && N0s) ? N0s->spec(ic+2,jc) : aEE;
    const float *oSS = (aSS && N0s) ? N0s->spec(ic,jc-2) : aSS, *oNN = (aNN && N0s) ? N0s->spec(ic,jc+2) : aNN;

    const double rdx = 1.0/p->DXN[IP];
    const double rdy = 1.0/p->DYN[JP];
    const double rdth = 1.0/g.dtheta;

    const double d = e->depth(ic,jc);
    const double U = e->U(ic,jc), V = e->V(ic,jc);
    const bool kin = (e->refr(ic,jc)==1);

    const double ddx = e->ddx(ic,jc), ddy = e->ddy(ic,jc);
    const double dUdx = e->dUdx(ic,jc), dUdy = e->dUdy(ic,jc), dVdx = e->dVdx(ic,jc), dVdy = e->dVdy(ic,jc);

    const int ma = m0[q], mb = m1[q], nq = mb-ma+1;

    if(src!=nullptr)
    src->compute(N,d,kc,cgc,P.data(),D.data());

    const double sig = g.sig[0];
    const double cgl = cgc[0];
    const double A = kin ? seastate_refraction(sig,kc[0],d) : 0.0;

    for(int n=0; n<nq; ++n)
    {
    const int m = ma+n;
    const int b = g.bin(0,m);
    const double cs = g.costh[m], sn = g.sinth[m];
    const int mm = (m-1+ndir)%ndir, mp = (m+1)%ndir;

    // c_theta at the faces m+1/2 and m-1/2
    double ctp=0.0, ctm=0.0;

        if(refraction && kin)
        {
        const double sp = g.sinthf[m],  cp_ = g.costhf[m];
        const double sm = g.sinthf[mm], cm  = g.costhf[mm];

        ctp = A*(sp*ddx - cp_*ddy) + sp*cp_*(dUdx-dVdy) + sp*sp*dVdx - cp_*cp_*dUdy;
        ctm = A*(sm*ddx - cm*ddy)  + sm*cm*(dUdx-dVdy)  + sm*sm*dVdx - cm*cm*dUdy;
        }

    // geographic face velocities
    const double cxc = cgl*cs + U;
    const double cyc = cgl*sn + V;

    const double cxw = W.cg  ? 0.5*(cxc + double(W.cg[0])*cs + W.U)   : cxc;
    const double cxe = E.cg  ? 0.5*(cxc + double(E.cg[0])*cs + E.U)   : cxc;
    const double cys = S.cg  ? 0.5*(cyc + double(S.cg[0])*sn + S.V)   : cyc;
    const double cyn = Nn.cg ? 0.5*(cyc + double(Nn.cg[0])*sn + Nn.V) : cyc;

    const double out  = (pos(cxe) - neg(cxw))*rdx + (pos(cyn) - neg(cys))*rdy;
    const double outt = (pos(ctp) - neg(ctm))*rdth;

    di[n] = rdt + th*(out + outt) + (src!=nullptr ? D[b] : 0.0);

    double r = rdt*double(Nc0[b]) - th1*(out + outt)*double(Nc0[b]);

    if(src!=nullptr)
    r += P[b];

    // zero-gradient sides: inflow of the cell's own spectrum
    double aself = 0.0;
    if(W.self)  aself += pos(cxw)*rdx;
    if(E.self)  aself -= neg(cxe)*rdx;
    if(S.self)  aself += pos(cys)*rdy;
    if(Nn.self) aself -= neg(cyn)*rdy;

        if(aself>0.0)
        {
        r += th1*aself*double(Nc0[b]);

            if(di[n]-th*aself>=0.5*di[n])
            di[n] -= th*aself;
            else
            r += th*aself*double(N[b]);
        }

    // inflow from the neighbours, new (latest) and old values
    if(W.N)  r += pos(cxw)*rdx*(th*W.N[b]  + th1*W0[b]);
    if(E.N)  r -= neg(cxe)*rdx*(th*E.N[b]  + th1*E0[b]);
    if(S.N)  r += pos(cys)*rdy*(th*S.N[b]  + th1*S0[b]);
    if(Nn.N) r -= neg(cyn)*rdy*(th*Nn.N[b] + th1*N0n[b]);

    // second-order correction of the geographic fluxes
    {
    const double fe = cxe>0.0 ? tvd(cxe,aW,N,aE,b)  : tvd(cxe,aEE,aE,N,b);
    const double fw = cxw>0.0 ? tvd(cxw,aWW,aW,N,b) : tvd(cxw,aE,N,aW,b);
    const double fn = cyn>0.0 ? tvd(cyn,aS,N,aN,b)  : tvd(cyn,aNN,aN,N,b);
    const double fs = cys>0.0 ? tvd(cys,aSS,aS,N,b) : tvd(cys,aN,N,aS,b);

    const double ge = cxe>0.0 ? tvd(cxe,oW,Nc0,oE,b)  : tvd(cxe,oEE,oE,Nc0,b);
    const double gw = cxw>0.0 ? tvd(cxw,oWW,oW,Nc0,b) : tvd(cxw,oE,Nc0,oW,b);
    const double gn = cyn>0.0 ? tvd(cyn,oS,Nc0,oN,b)  : tvd(cyn,oNN,oN,Nc0,b);
    const double gs = cys>0.0 ? tvd(cys,oSS,oS,Nc0,b) : tvd(cys,oN,Nc0,oS,b);

    r -= th*((fe - fw)*rdx + (fn - fs)*rdy) + th1*((ge - gw)*rdx + (gn - gs)*rdy);
    }

    // theta: within the quadrant implicit (tridiagonal), outside with the latest values; old values explicit
    r += th1*(pos(ctm)*double(Nc0[g.bin(0,mm)]) - neg(ctp)*double(Nc0[g.bin(0,mp)]))*rdth;

    la[n] = up[n] = 0.0;

        if(n>0)
        la[n] = -th*pos(ctm)*rdth;
        else
        r += th*pos(ctm)*rdth*N[g.bin(0,mm)];

        if(n<nq-1)
        up[n] = th*neg(ctp)*rdth;
        else
        r -= th*neg(ctp)*rdth*N[g.bin(0,mp)];

    rhs[n] = r;
    }

    // Thomas algorithm
    for(int n=0; n<nq; ++n)
    {
    const double den = di[n] - (n>0 ? la[n]*cp[n-1] : 0.0);

        if(!(den>1.0e-300))
        {
        for(int nn=0; nn<nq; ++nn)
        N[g.bin(0,ma+nn)] = 0.0f;

        i = ic;
        j = jc;
        return;
        }

    cp[n] = up[n]/den;
    dp[n] = (rhs[n] - (n>0 ? la[n]*dp[n-1] : 0.0))/den;
    }

    sol[nq-1] = dp[nq-1];
    for(int n=nq-2; n>=0; --n)
    sol[n] = dp[n] - cp[n]*sol[n+1];

    for(int n=0; n<nq; ++n)
    N[g.bin(0,ma+n)] = float(std::max(sol[n],0.0));

    i = ic;
    j = jc;
}
