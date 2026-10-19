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

//  REEFAMR: the patches of the mesh refinement in checkpoints (P 48, I 45), for every module with its
//  list of arrays per patch (fields(patch, visit): visit(name, kind, double*, int*), kind S slice,
//  F node field FIJK, C cell field IJK of the patch lexer; a nullptr array is skipped).  The blocks
//  1.. of a multi-block grid (level-0 patches, each whole on one rank) are patches here as well.
//
//  Writing: every rank its patches as a rank blob <tag>: per patch its level and box, the ranges of
//  its lexer and the arrays (ghost and EXT cells included).
//  Reading: ini has built the hierarchy for this decomposition (its patches can be cut differently
//  from the checkpoint's).  A patch with the same level and box as one of the checkpoint takes its
//  arrays as they are (same decomposition: bitwise); otherwise every column (I,J) of the patch, ghost
//  and EXT cells included, takes its values from the checkpoint patch of the level whose box holds
//  (I,J), else from one whose lexer reaches it.  The refinement must be the same (static boxes, or
//  the layout of the checkpoint).

#include"reefamr.h"
#include"reefckpt.h"
#include"lexer.h"
#include<vector>
#include<string>

using namespace std;

namespace
{
struct ckpatch
{
    int lev, I0, I1, J0, J1, imin, imax, jmin, jmax, kmaxF, kmax;
    vector<const double*> rec;      // per array of the list, into the blob
    vector<long long> n;
};

long long kcount(char k, int kmaxF, int kmax)
{
    return k=='S' ? 1 : (k=='F' ? kmaxF : kmax);
}
}

void reefamr::ckpt_write_patches(reefckpt &ck, const char *tag, const ckpt_fields_fn &fields)
{
    vector<double> b;

    for(auto q : P)
    {
        if(q->lev<1 && q->blk==0)
        continue;

        lexer *pp = q->pp;

        int nrec=0;
        fields(q,[&](const char*,char,double*,int*){++nrec;});

        const double hd[12] = {double(q->lev),double(q->I0),double(q->I1),double(q->J0),double(q->J1),
                               double(pp->imin),double(pp->imax),double(pp->jmin),double(pp->jmax),double(pp->kmaxF),double(pp->kmax),double(nrec)};
        b.insert(b.end(),hd,hd+12);

        const long long ncol = (long long)pp->imax*pp->jmax;
        auto put = [&](char k, auto *V)
        {
            const long long n = ncol*kcount(k,pp->kmaxF,pp->kmax);
            b.push_back(V==nullptr ? -1.0 : double(n));
            if(V!=nullptr)
            for(long long m=0; m<n; ++m)
            b.push_back(double(V[m]));
        };
        fields(q,[&](const char*, char k, double *V, int *I){ if(V!=nullptr) put(k,V); else put(k,I); });
    }

    ck.rank_blob(tag,b);
}

bool reefamr::ckpt_read_patches(reefckpt &ck, const char *tag, const ckpt_fields_fn &fields, string &msg)
{
    // the patches of all pieces of the checkpoint
    vector<vector<double>> blobs = ck.get_rank_blobs(tag);
    vector<ckpatch> cp;

    for(const auto &bl : blobs)
    {
        size_t q=0;
        while(q+12<=bl.size())
        {
            ckpatch c;
            c.lev = int(bl[q]); c.I0 = int(bl[q+1]); c.I1 = int(bl[q+2]); c.J0 = int(bl[q+3]); c.J1 = int(bl[q+4]);
            c.imin = int(bl[q+5]); c.imax = int(bl[q+6]); c.jmin = int(bl[q+7]); c.jmax = int(bl[q+8]);
            c.kmaxF = int(bl[q+9]); c.kmax = int(bl[q+10]);
            const int nrec = int(bl[q+11]);
            q += 12;

            // the arrays in the order of fnpf_ckpt_fields
            for(int r=0; r<nrec; ++r)
            {
                if(q>=bl.size())
                {
                    msg = "the mesh refinement in the checkpoint is incomplete";
                    return false;
                }
                const long long n = (long long)bl[q++];
                c.n.push_back(n);
                c.rec.push_back(n>=0 ? &bl[q] : nullptr);
                if(n>0)
                q += size_t(n);
                if(q>bl.size())
                {
                    msg = "the mesh refinement in the checkpoint is incomplete";
                    return false;
                }
            }
            cp.push_back(c);
        }
    }

    // the patches of this rank
    for(auto q : P)
    {
        if(q->lev<1 && q->blk==0)
        continue;

        lexer *pp = q->pp;
        const long long ncol = (long long)pp->imax*pp->jmax;

        // the arrays of the patch, in the order of fnpf_ckpt_fields
        vector<pair<char,void*>> arr;
        vector<char> isint;
        fields(q,[&](const char*, char k, double *V, int *I)
        {
            if(V!=nullptr || I==nullptr) { arr.push_back({k,(void*)V}); isint.push_back(0); }
            else { arr.push_back({k,(void*)I}); isint.push_back(1); }
        });

        // the same patch in the checkpoint
        int same = -1;
        for(int s=0; s<(int)cp.size() && same<0; ++s)
        {
            const ckpatch &c = cp[s];
            if(c.lev==q->lev && c.I0==q->I0 && c.I1==q->I1 && c.J0==q->J0 && c.J1==q->J1
            && c.imax==pp->imax && c.jmax==pp->jmax && c.kmaxF==pp->kmaxF && c.kmax==pp->kmax)
            same = s;
        }

        // per column the source patch and its column, else -1
        vector<int> src((size_t)ncol,-1);
        vector<long long> scol((size_t)ncol,0);
        if(same>=0)
        {
            for(long long m=0; m<ncol; ++m)
            {
                src[m] = same;
                scol[m] = m;
            }
        }
        else
        for(int ii=pp->imin; ii<pp->imin+pp->imax; ++ii)
        for(int jj=pp->jmin; jj<pp->jmin+pp->jmax; ++jj)
        {
            const int I = ii + q->I0 - EXT, J = jj + q->J0 - EXT;
            int best = -1;
            for(int s=0; s<(int)cp.size(); ++s)
            {
                const ckpatch &c = cp[s];
                if(c.lev!=q->lev || c.kmaxF!=pp->kmaxF || c.kmax!=pp->kmax)
                continue;
                const int si = I - c.I0 + EXT, sj = J - c.J0 + EXT;
                if(si<c.imin || si>=c.imin+c.imax || sj<c.jmin || sj>=c.jmin+c.jmax)
                continue;
                const bool inside = (I>=c.I0 && I<=c.I1 && J>=c.J0 && J<=c.J1);
                if(inside || best<0)
                best = s;
                if(inside)
                break;
            }
            if(best>=0)
            {
                const ckpatch &c = cp[best];
                const long long m = (long long)(ii-pp->imin)*pp->jmax + (jj-pp->jmin);
                src[m] = best;
                scol[m] = (long long)((I - c.I0 + EXT) - c.imin)*c.jmax + ((J - c.J0 + EXT) - c.jmin);
            }
        }

        // copy
        for(size_t r=0; r<arr.size(); ++r)
        {
            if(arr[r].second==nullptr)
            continue;
            const long long nk = kcount(arr[r].first,pp->kmaxF,pp->kmax);
            for(long long m=0; m<ncol; ++m)
            {
                const int s = src[m];
                if(s<0 || r>=cp[s].rec.size() || cp[s].rec[r]==nullptr)
                continue;
                const double *from = cp[s].rec[r] + scol[m]*nk;
                if(isint[r])
                {
                    int *to = (int*)arr[r].second + m*nk;
                    for(long long kk=0; kk<nk; ++kk)
                    to[kk] = int(from[kk]);
                }
                else
                {
                    double *to = (double*)arr[r].second + m*nk;
                    for(long long kk=0; kk<nk; ++kk)
                    to[kk] = from[kk];
                }
            }
        }
    }

    return true;
}

// the layout of the hierarchy: the tile maps and the hysteresis ages of every level, the global patch
// table (ranks of the placement), the regrid counters
void reefamr::ckpt_write_layout(reefckpt &ck)
{
    vector<double> b;
    b.push_back(maxlev);
    b.push_back(p0->mpi_size);
    b.push_back(regrids);
    b.push_back(regrids_skipped);
    for(int l=1; l<=maxlev; ++l)
    {
        const size_t n = (l<(int)gtile.size()) ? gtile[l].size() : 0;
        b.push_back(double(n));
        for(size_t t=0; t<n; ++t)
        b.push_back(gtile[l][t]);
        const size_t na = (l<(int)tage.size()) ? tage[l].size() : 0;
        b.push_back(double(na));
        for(size_t t=0; t<na; ++t)
        b.push_back(tage[l][t]);
    }
    b.push_back(double(GP.size()));
    for(const auto &g : GP)
    {
        b.push_back(g.lev); b.push_back(g.I0); b.push_back(g.I1); b.push_back(g.J0); b.push_back(g.J1);
        b.push_back(g.rank); b.push_back(g.lid);
    }
    ck.blob("amr_layout",b);
}

bool reefamr::ckpt_read_layout(reefckpt &ck, string &msg)
{
    const vector<double> b = ck.get_blob("amr_layout");
    size_t q=0;
    auto nx = [&]() { return q<b.size() ? b[q++] : -1.0; };

    if(b.empty() || int(nx())!=maxlev)
    {
        msg = "the mesh refinement of the checkpoint has other levels (G 1)";
        return false;
    }
    const int nranks = int(nx());
    regrids = int(nx());
    regrids_skipped = int(nx());

    layoutM.assign(maxlev+1,vector<unsigned char>());
    if((int)tage.size()<maxlev+1)
    tage.resize(maxlev+1);
    for(int l=1; l<=maxlev; ++l)
    {
        const size_t n = size_t(nx());
        layoutM[l].resize(n);
        for(size_t t=0; t<n; ++t)
        layoutM[l][t] = (unsigned char)nx();
        const size_t na = size_t(nx());
        tage[l].resize(na);
        for(size_t t=0; t<na; ++t)
        tage[l][t] = (unsigned char)nx();
    }

    // the ranks of the patches: kept on the same number of ranks (same decomposition)
    const size_t ng = size_t(nx());
    vector<reefamr_gpatch> gp(ng);
    for(auto &g : gp)
    {
        g.lev = int(nx()); g.I0 = int(nx()); g.I1 = int(nx()); g.J0 = int(nx()); g.J1 = int(nx());
        g.rank = int(nx()); g.lid = int(nx());
    }
    if(q>b.size())
    {
        msg = "the layout of the mesh refinement in the checkpoint is incomplete";
        return false;
    }
    if(nranks==p0->mpi_size && ck.same_decomposition())
    {
        GPckpt = gp;
        keep_ranks = true;
    }
    return true;
}
