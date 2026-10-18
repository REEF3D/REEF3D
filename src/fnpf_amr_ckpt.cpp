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

//  FNPF mesh refinement in checkpoints (P 48, I 45).
//
//  Writing: every rank its patches as a rank blob "fnpf_amr": per patch its level and box, the ranges
//  of its lexer and the state of fnpf_ckpt_fields (ghost and EXT cells included).
//  Reading: ini has built the hierarchy for this decomposition (its patches can be cut differently
//  from the checkpoint's).  A patch with the same level and box as one of the checkpoint takes its
//  arrays as they are (same decomposition: bitwise); otherwise every column (I,J) of the patch, ghost
//  and EXT cells included, takes its values from the checkpoint patch of the level whose box holds
//  (I,J), else from one whose lexer reaches it.  The refinement must be the same (static boxes, or the
//  same layout); the regrids of moving bodies (G 12) need the 6DOF state, not there yet.

#include"fnpf_amr.h"
#include"reefckpt.h"
#include"fnpf_ckpt_fields.h"
#include<vector>
#include<string>

using namespace std;

namespace
{
struct ckpatch
{
    int lev, I0, I1, J0, J1, imin, imax, jmin, jmax, kmaxF, kmax;
    vector<const double*> rec;      // per array of fnpf_ckpt_fields, into the blob
    vector<long long> n;
};

long long kcount(char k, int kmaxF, int kmax)
{
    return k=='S' ? 1 : (k=='F' ? kmaxF : kmax);
}
}

void fnpf_amr::ckpt_write(reefckpt &ck)
{
    vector<double> b;

    for(auto q : P)
    {
        if(q->lev<1)
        continue;

        fnpf_amr_patch *pc = FP(q);
        lexer *pp = pc->pp;

        int nrec=0;
        fnpf_ckpt_fields(pp,pc->c,[&](const char*,char,double*){++nrec;},[&](const char*,char,int*){++nrec;});

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
        fnpf_ckpt_fields(pp,pc->c,
            [&](const char*, char k, double *V){ put(k,V); },
            [&](const char*, char k, int *V){ put(k,V); });
    }

    ck.rank_blob("fnpf_amr",b);
}

bool fnpf_amr::ckpt_read(reefckpt &ck, string &msg)
{
    // the patches of all pieces of the checkpoint
    vector<vector<double>> blobs = ck.get_rank_blobs("fnpf_amr");
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
        if(q->lev<1)
        continue;

        fnpf_amr_patch *pc = FP(q);
        lexer *pp = pc->pp;
        const long long ncol = (long long)pp->imax*pp->jmax;

        // the arrays of the patch, in the order of fnpf_ckpt_fields
        vector<pair<char,void*>> arr;
        vector<char> isint;
        fnpf_ckpt_fields(pp,pc->c,
            [&](const char*, char k, double *V){ arr.push_back({k,(void*)V}); isint.push_back(0); },
            [&](const char*, char k, int *V){ arr.push_back({k,(void*)V}); isint.push_back(1); });

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
