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

//  CFD mesh refinement in checkpoints (P 48, I 45).
//
//  The hierarchy: the tile maps of the levels, the patch table (boxes, ranks), the hysteresis of the
//  tiles and the counters of the refinement (blob cfd_amr_layout).  A restart regrids to these tile
//  maps (ck_seed; on the same decomposition with the ranks of the checkpoint, GPforce).
//  The patches: every array of a patch (its fields, slices, solver arrays: ckpt_registry.h with the
//  lexer of the patch) with its box (rank blob cfd_amr_patches).  A patch of the same box takes them
//  as they are (bitwise); another one (a changed decomposition, G 40 units cut anew) takes its cells
//  from the patches of its level in the checkpoint that hold them, the solver arrays start again.

#include"cfd_amr.h"
#include"lexer.h"
#include"reefckpt.h"
#include"ckpt_registry.h"
#include"amr_load.h"
#include<cstring>
#include<fstream>
#include<iomanip>
#include<sys/stat.h>

using namespace std;

namespace
{
// index of a cell in the arrays of lexer q (field layout)
inline int cix(const lexer *q, int i, int j, int k)
{
    return (i-q->imin)*q->jmax*q->kmax + (j-q->jmin)*q->kmax + k-q->kmin;
}

struct ckpatch
{
    int shadow, lev, lo[3], hi[3];
    int imin, imax, jmin, jmax, kmin, kmax;
    double dt_l, dt, dt_old, simtime;    // the step of the level, the time values of the patch lexer
    vector<int> code;
    vector<vector<double>> data;
};

double code_of(const ckpt_registry::entry &e)
{
    return (e.kind=='C' ? 0.0 : (e.kind=='S' ? 2.0 : 8.0)) + (e.isint ? 1.0 : 0.0) + (e.local ? 4.0 : 0.0);
}

long long count_of(const ckpt_registry::entry &e, const lexer *q)
{
    if(e.kind=='R')
    return e.n;
    const long long s = (long long)q->imax*q->jmax;
    return e.kind=='C' ? s*q->kmax : s;
}
}

void cfd_amr::ckpt_write(reefckpt &ck)
{
    // ---- the hierarchy (the same on every rank)
    vector<double> L;
    L.push_back(maxlev);
    L.push_back(nranks);
    for(int l=1; l<=maxlev; ++l)
    {
        const size_t n = (l<(int)tmask.size()) ? tmask[l].size() : 0;
        L.push_back(double(n));
        for(size_t t=0; t<n; ++t)
        L.push_back(tmask[l][t]);
    }
    L.push_back(double(GP.size()));
    for(const auto &g : GP)
    {
        L.push_back(g.lev);
        for(int d=0; d<3; ++d) L.push_back(g.lo[d]);
        for(int d=0; d<3; ++d) L.push_back(g.hi[d]);
        L.push_back(g.rank); L.push_back(g.lid); L.push_back(g.unit);
    }
    L.push_back(double(tage.size()));
    for(const auto &t : tage)
    {
        L.push_back(double(t.size()));
        for(int v : t)
        L.push_back(v);
    }
    L.push_back(double(nregrid));
    L.push_back(vol_step); L.push_back(vol_net); L.push_back(vol_ini ? 1.0 : 0.0);
    L.push_back(double(th_end.size()));
    for(double v : th_end)
    L.push_back(v);
    L.push_back(wpatch);
    L.push_back(press_step ? 1.0 : 0.0);
    ck.blob("cfd_amr_layout",L);

    // ---- the arrays of the patches of this rank
    vector<double> B;
    for(r3patch *c : P)
    {
        const lexer *q = c->pp;
        vector<ckpt_registry::entry> e;
        ckpt_registry::visit(q,[&](const ckpt_registry::entry &x){ e.push_back(x); });

        B.push_back(c->shadow ? 1.0 : 0.0);
        B.push_back(c->lev);
        for(int d=0; d<3; ++d) B.push_back(c->lo[d]);
        for(int d=0; d<3; ++d) B.push_back(c->hi[d]);
        B.push_back(q->imin); B.push_back(q->imax); B.push_back(q->jmin); B.push_back(q->jmax); B.push_back(q->kmin); B.push_back(q->kmax);
        B.push_back(c->shadow ? 0.0 : static_cast<cfd_amr_patch*>(c)->dt_l);
        B.push_back(q->dt); B.push_back(q->dt_old); B.push_back(q->simtime);
        B.push_back(double(e.size()));
        for(const auto &x : e)
        {
            const long long n = count_of(x,q);
            B.push_back(code_of(x));
            B.push_back(double(n));
            if(x.isint)
            for(long long m=0; m<n; ++m)
            B.push_back(double(((const int*)x.V)[m]));
            else
            B.insert(B.end(),(const double*)x.V,(const double*)x.V+n);
        }
    }
    ck.rank_blob("cfd_amr_patches",B);

    // G 7 1: the start of the synchronisation projection (the correction of the last one, per grid)
    vector<double> S;
    S.push_back(sync_x_ok ? 1.0 : 0.0);
    S.push_back(sync_bref);
    S.push_back(double(sync_x.size()));
    for(const auto &x : sync_x)
    {
        S.push_back(double(x.size()));
        S.insert(S.end(),x.begin(),x.end());
    }
    ck.rank_blob("cfd_amr_sync",S);
}

bool cfd_amr::ckpt_read(reefckpt &ck, string &msg)
{
    // ---- the hierarchy
    const vector<double> L = ck.get_blob("cfd_amr_layout");
    size_t k=0;
    auto nx = [&]() { return k<L.size() ? L[k++] : -1.0; };

    if(L.empty() || int(nx())!=maxlev)
    {
        msg = "the mesh refinement of the checkpoint has other levels (G 1)";
        return false;
    }
    const int ckranks = int(nx());

    vector<vector<unsigned char>> M(maxlev+1);
    for(int l=1; l<=maxlev; ++l)
    {
        const size_t n = size_t(nx());
        M[l].resize(n);
        for(size_t t=0; t<n; ++t)
        M[l][t] = (unsigned char)nx();
        if(l<(int)tmask.size() && n!=tmask[l].size())
        {
            msg = "the tiles of the mesh refinement in the checkpoint differ (G 4)";
            return false;
        }
    }
    vector<r3gpatch> G((size_t)nx());
    for(auto &g : G)
    {
        g.lev = int(nx());
        for(int d=0; d<3; ++d) g.lo[d] = int(nx());
        for(int d=0; d<3; ++d) g.hi[d] = int(nx());
        g.rank = int(nx()); g.lid = int(nx()); g.unit = int(nx());
    }
    vector<vector<int>> T((size_t)nx());
    for(auto &t : T)
    {
        t.resize((size_t)nx());
        for(auto &v : t)
        v = int(nx());
    }
    const long ck_nregrid = long(nx());
    const double ck_vs = nx(), ck_vn = nx();
    const bool ck_vi = nx()>0.5;
    vector<double> th((size_t)nx());
    for(auto &v : th)
    v = nx();
    const double ck_wpatch = nx();
    const double ck_pstep = nx();
    if(k>L.size())
    {
        msg = "the mesh refinement in the checkpoint is incomplete";
        return false;
    }

    // ---- regrid to the tile maps of the checkpoint (the ranks of its patches on the same decomposition)
    const bool same = ck.same_decomposition() && ckranks==nranks;
    ck_seed = M;
    GPforce = same ? &G : nullptr;
    regrid(p0,pgc0,false);
    ck_seed.clear();
    GPforce = nullptr;

    tage = T;
    nregrid = ck_nregrid;
    vol_step = ck_vs; vol_net = ck_vn; vol_ini = ck_vi;
    if(th.size()==th_end.size())
    th_end = th;
    if(ck_wpatch>0.0)
    wpatch = ck_wpatch;
    if(ck_pstep>=0.0)
    press_step = ck_pstep>0.5;

    // ---- the arrays of the patches of all pieces
    vector<ckpatch> cp;
    for(const auto &bl : ck.get_rank_blobs("cfd_amr_patches"))
    {
        size_t q=0;
        while(q+21<=bl.size())
        {
            ckpatch c;
            c.shadow = int(bl[q++]); c.lev = int(bl[q++]);
            for(int d=0; d<3; ++d) c.lo[d] = int(bl[q++]);
            for(int d=0; d<3; ++d) c.hi[d] = int(bl[q++]);
            c.imin = int(bl[q++]); c.imax = int(bl[q++]); c.jmin = int(bl[q++]); c.jmax = int(bl[q++]); c.kmin = int(bl[q++]); c.kmax = int(bl[q++]);
            c.dt_l = bl[q++];
            c.dt = bl[q++]; c.dt_old = bl[q++]; c.simtime = bl[q++];
            const int ne = int(bl[q++]);
            for(int r=0; r<ne; ++r)
            {
                if(q+2>bl.size())
                {
                    msg = "the mesh refinement in the checkpoint is incomplete";
                    return false;
                }
                c.code.push_back(int(bl[q++]));
                const size_t n = size_t(bl[q++]);
                if(q+n>bl.size())
                {
                    msg = "the mesh refinement in the checkpoint is incomplete";
                    return false;
                }
                c.data.emplace_back(bl.begin()+q,bl.begin()+q+n);
                q += n;
            }
            cp.push_back(std::move(c));
        }
    }

    int exact=0, mapped=0;
    for(r3patch *c : P)
    {
        lexer *q = c->pp;
        vector<ckpt_registry::entry> e;
        ckpt_registry::visit(q,[&](const ckpt_registry::entry &x){ e.push_back(x); });

        // a patch of the same box: every array as it is
        const ckpatch *m = nullptr;
        for(const auto &o : cp)
        if(o.shadow==(c->shadow?1:0) && o.lev==c->lev && o.lo[0]==c->lo[0] && o.lo[1]==c->lo[1] && o.lo[2]==c->lo[2]
           && o.hi[0]==c->hi[0] && o.hi[1]==c->hi[1] && o.hi[2]==c->hi[2] && o.code.size()==e.size())
        {
            bool fit = true;
            for(size_t r=0; r<e.size() && fit; ++r)
            if(o.code[r]!=int(code_of(e[r])) || (long long)o.data[r].size()!=count_of(e[r],q))
            fit = false;
            if(fit)
            {
                m = &o;
                break;
            }
        }

        // the step of the patch (the previous step of its level: dt_old of its first step)
        if(!c->shadow)
        for(const auto &o : cp)
        if(!o.shadow && o.lev==c->lev)
        {
            const ckpatch &src = (m!=nullptr) ? *m : o;
            static_cast<cfd_amr_patch*>(c)->dt_l = src.dt_l;
            q->dt = src.dt;
            q->dt_old = src.dt_old;
            q->simtime = src.simtime;
            break;
        }

        if(m!=nullptr)
        {
            for(size_t r=0; r<e.size(); ++r)
            {
                if(e[r].local && !same)
                continue;
                const vector<double> &v = m->data[r];
                if(e[r].isint)
                for(size_t t=0; t<v.size(); ++t)
                ((int*)e[r].V)[t] = int(v[t]);
                else
                std::copy(v.begin(),v.end(),(double*)e[r].V);
            }
            ++exact;
            continue;
        }

        if(c->shadow)
        continue;

        // another box: the cells (I,J,K) of the level from the patches of the checkpoint that hold them,
        // the interior of a patch first
        const int jk = q->jmax*q->kmax;
        for(size_t r=0; r<e.size(); ++r)
        {
            // the scalars of the objects (the same on all patches of the level): from one of them
            if(e[r].kind=='R' && !e[r].local)
            {
                for(const auto &o : cp)
                if(!o.shadow && o.lev==c->lev && o.code.size()==e.size() && o.code[r]==int(code_of(e[r])) && (long long)o.data[r].size()==e[r].n)
                {
                    std::copy(o.data[r].begin(),o.data[r].end(),(double*)e[r].V);
                    break;
                }
                continue;
            }
            if(e[r].kind=='R' || e[r].local)
            continue;
            const bool cell = (e[r].kind=='C');

            for(int ii=q->imin; ii<q->imin+q->imax; ++ii)
            for(int jj=q->jmin; jj<q->jmin+q->jmax; ++jj)
            for(int kk=(cell ? q->kmin : 0); kk<(cell ? q->kmin+q->kmax : 1); ++kk)
            {
                const int I = c->lo[0]+ii, J = c->lo[1]+jj, K = c->lo[2]+kk;

                const ckpatch *src = nullptr;
                int best = 2;
                for(const auto &o : cp)
                {
                    if(o.shadow || o.lev!=c->lev || o.code.size()!=e.size() || o.code[r]!=int(code_of(e[r])))
                    continue;
                    const int oi = I-o.lo[0], oj = J-o.lo[1], ok = K-o.lo[2];
                    if(oi<o.imin || oi>=o.imin+o.imax || oj<o.jmin || oj>=o.jmin+o.jmax)
                    continue;
                    if(cell && (ok<o.kmin || ok>=o.kmin+o.kmax))
                    continue;
                    const bool inner = (I>=o.lo[0] && I<=o.hi[0] && J>=o.lo[1] && J<=o.hi[1] && (!cell || (K>=o.lo[2] && K<=o.hi[2])));
                    const int rank = inner ? 0 : 1;
                    if(rank<best)
                    {
                        best = rank;
                        src = &o;
                        if(rank==0)
                        break;
                    }
                }
                if(src==nullptr)
                continue;

                const int oi = I-src->lo[0], oj = J-src->lo[1], ok = K-src->lo[2];
                const size_t si = cell ? (size_t)(oi-src->imin)*src->jmax*src->kmax + (size_t)(oj-src->jmin)*src->kmax + (ok-src->kmin)
                                       : (size_t)(oi-src->imin)*src->jmax + (oj-src->jmin);
                const size_t di = cell ? (size_t)(ii-q->imin)*jk + (size_t)(jj-q->jmin)*q->kmax + (kk-q->kmin)
                                       : (size_t)(ii-q->imin)*q->jmax + (jj-q->jmin);
                const double v = src->data[r][si];
                if(e[r].isint)
                ((int*)e[r].V)[di] = int(v);
                else
                ((double*)e[r].V)[di] = v;
            }
        }
        ++mapped;
    }

    // G 7 1: the start of the synchronisation projection (the same decomposition: the same grids)
    if(same)
    {
        const vector<double> S = ck.get_rank_blob("cfd_amr_sync",p0->mpirank);
        if(S.size()>=3)
        {
            size_t q=0;
            const bool ok = S[q++]>0.5;
            const double bref = S[q++];
            vector<vector<double>> X(size_t(S[q++]));
            bool fit = true;
            for(auto &x : X)
            {
                const size_t n = (q<S.size()) ? size_t(S[q++]) : 0;
                if(q+n>S.size()) { fit = false; break; }
                x.assign(S.begin()+q,S.begin()+q+n);
                q += n;
            }
            if(fit)
            {
                sync_x = X;
                sync_x_ok = ok;
                sync_bref = bref;
            }
        }
    }

    if(p0->mpirank==0 && verbose)
    cout<<"CFD AMR: restart from the checkpoint, "<<exact<<" patches as they were, "<<mapped<<" mapped (rank 0)"<<endl;

    return true;
}

// the load map (G 45, G 46): per step of the window of G 46 the patch cells over every level-0 cell
// (G 7 1: times the steps of their level), from the patch table of all ranks (rank 0), and the level-0
// cells (not solid) and patch cells of this rank
void cfd_amr::lmap_step(lexer *p)
{
    if(p->G45<=0 && p->G46<=0.0)
    return;

    // a new window of G 46 (the rebalance checks every G 47 steps)
    const int G47 = MAX(p->G47,1);
    const int win = (p->count-1)/G47;
    if(win!=lm_win)
    {
        lm_win = win;
        lm_acc.clear();
        lm_steps = 0;
        lm_n0 = lm_np = 0.0;
    }

    long act=0;
    for(int i=0; i<p0->knox; ++i)
    for(int j=0; j<p0->knoy; ++j)
    for(int k=0; k<p0->knoz; ++k)
    if(p0->flag4[cix(p0,i,j,k)]>0)
    ++act;
    lm_n0 += double(act);
    for(int qi=0; qi<npat; ++qi)
    lm_np += ((p0->G7==1) ? double(1<<P[qi]->lev) : 1.0)*double(P[qi]->n[0])*P[qi]->n[1]*P[qi]->n[2];
    ++lm_steps;

    const long long NX = gn0[0], NY = gn0[1], NZ = gn0[2];
    if(p->mpirank==0)
    {
    if(lm_nregrid!=nregrid || lm_cur.size()!=(size_t)(NX*NY*NZ))
    {
        lm_nregrid = nregrid;
        lm_cur.assign((size_t)(NX*NY*NZ),0.0);
        for(const auto &g : GP)
        {
            const double w = (p0->G7==1) ? double(1<<g.lev) : 1.0;
            int R[3];
            for(int d=0; d<3; ++d)
            R[d] = ratio(g.lev,d);
            for(int I=g.lo[0]; I<=g.hi[0]; ++I)
            for(int J=g.lo[1]; J<=g.hi[1]; ++J)
            for(int K=g.lo[2]; K<=g.hi[2]; ++K)
            {
                const long long a = I/R[0], b = J/R[1], c = K/R[2];
                if(a>=0 && a<NX && b>=0 && b<NY && c>=0 && c<NZ)
                lm_cur[(size_t)((a*NY+b)*NZ+c)] += w;
            }
        }
    }
    if(lm_acc.size()!=lm_cur.size())
    lm_acc.assign(lm_cur.size(),0.0);
    for(size_t n=0; n<lm_cur.size(); ++n)
    lm_acc[n] += lm_cur[n];
    }

    // G 45: every G 45 steps (collective)
    if(p->G45>0 && p->count%p->G45==0)
    write_loadmap(p);
}

// the load map (G 45, G 46): the work per step of every level-0 cell, 1 plus the patch cells over it
// (averaged over the window of G 46) times the weight of a patch cell.  The weight: fitted to the
// busy times of the ranks in the window (busy = a (level-0 cells + w patch cells), least squares over
// the ranks), else wpatch (G 41, or measured at the regrids).  Rank 0 writes it for the weighted
// partition (DIVEMesh M 20 5, M 35 1; REEF3D M 41 1), ./REEF3D_CFD_AMR/REEF3D_CFD_AMR_loadmap.dat:
// nx ny nz, then the cells (i*ny + j)*nz + k
void cfd_amr::write_loadmap(lexer *p)
{
    // the weight of a patch cell from the busy times of the ranks (collective)
    double w = (p->G41>0.0) ? p->G41 : wpatch;
    bool fit = false;
    if(p->G41<=0.0 && nranks>1 && lm_steps>0)
    {
        const double t = amr_window_busy(), x = lm_n0/double(lm_steps), y = lm_np/double(lm_steps);
        double v[6] = {x*x, x*y, y*y, t*x, t*y, t>0.0 ? 1.0 : 0.0}, s[6];
        MPI_Allreduce(v,s,6,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        const double det = s[0]*s[2]-s[1]*s[1];
        if(s[5]>=nranks && det>1.0e-6*s[0]*s[2])
        {
            const double A = (s[3]*s[2]-s[4]*s[1])/det, B = (s[0]*s[4]-s[1]*s[3])/det;
            if(A>0.0 && B>0.0)
            {
                w = MIN(MAX(B/A,0.25),16.0);
                fit = true;
            }
        }
    }

    if(p->mpirank!=0)
    return;

    const long long NX = gn0[0], NY = gn0[1], NZ = gn0[2];
    const double st = double(MAX(lm_steps,1L));
    vector<double> W((size_t)(NX*NY*NZ),1.0);
    if(lm_acc.size()==W.size())
    for(size_t n=0; n<W.size(); ++n)
    W[n] += w*lm_acc[n]/st;

    mkdir("./REEF3D_CFD_AMR",0777);
    ofstream out("./REEF3D_CFD_AMR/REEF3D_CFD_AMR_loadmap.dat");
    out<<"# REEF3D CFD AMR load map (G 45/46): work per step of every level-0 cell (1 + its patch cells times "<<w
       <<(fit ? ", measured" : "")<<"), averaged over "<<lm_steps<<" steps, count "<<p->count<<endl;
    out<<NX<<" "<<NY<<" "<<NZ<<endl;
    out<<setprecision(8);
    size_t n=0;
    for(double v : W)
    {
        out<<v<<((++n)%10==0 ? "\n" : " ");
    }
    out<<endl;

    const streamsize pr = cout.precision();
    cout<<"CFD AMR: load map, a patch cell "<<setprecision(3)<<w<<" level-0 cells"<<(fit ? " (fitted to the busy times of the ranks)" : "")<<endl;
    cout.precision(pr);
}

void cfd_amr::flush_outputs()
{
}
