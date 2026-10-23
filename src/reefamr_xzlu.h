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

#ifndef REEFAMR_XZLU_H_
#define REEFAMR_XZLU_H_

#include"reefmg_core.h"
#include<vector>
#include<algorithm>
#include<mpi.h>

//  Exact solve of the 7-point operator of a reefmg_core fine level with one cell in y (a patch of a
//  2D grid: an x-z problem), in place of its V-cycle.  The unknowns are ordered column by column
//  (k fastest), so the matrix is banded with nz sub- and super-diagonals: LU without pivoting in
//  band storage (the rows are those the line smoother of the V-cycle factorises without pivoting
//  too), O(nx nz^2) per factorisation, O(nx nz) per solve.  The couplings beyond the level (the
//  patch edge, below the bed node, above the top node) are 0 in the rows or left out, rows with
//  act 0 are identity rows.  factor_mpi: level 0 of a 2D grid partitioned in x, the problem of all
//  ranks factorised on every rank, the right-hand side gathered (one MPI_Allgatherv per solve).
struct reefamr_xzlu
{
    int nx = 0, nz = 0, w = 0;         // nx: of the whole problem (all ranks); w = 2 nz + 1, the band width
    std::vector<double> a;              // row m: entries of the columns m-nz .. m+nz

    // several ranks (factor_mpi): the levels of the ranks side by side in x, every rank factorises
    // and solves the whole problem (it is small: a 2D grid), the right-hand side gathered
    bool dist = false;
    MPI_Comm comm = MPI_COMM_NULL;
    int nxl = 0, i0 = 0;                // this rank: its columns, the first of them
    std::vector<int> cnt, dsp;          // per rank: unknowns and offset

    bool usable(const sc_level &L) const { return L.ny==1; }

    void factor(const sc_level &L, int prec)
    {
        dist = false;
        nxl = L.nx; i0 = 0;
        std::vector<double> G;
        pack(L,prec,G);
        factor_from(G,L.nx,L.nz);
    }

    // a 2D grid partitioned in x over the ranks of comm (this rank: the columns gi0 .. gi0+L.nx-1 of
    // gnx).  Collective; false on all ranks (nothing factorised) unless the ranks tile 0..gnx-1 with
    // one cell in y and the same nz
    bool factor_mpi(const sc_level &L, int prec, MPI_Comm cm, int gi0, int gnx)
    {
        int np=1;
        MPI_Comm_size(cm,&np);
        int me[3] = {gi0, L.nx, (L.ny==1) ? L.nz : -1};
        std::vector<int> all(3*np);
        MPI_Allgather(me,3,MPI_INT,all.data(),3,MPI_INT,cm);

        std::vector<int> cover(gnx,0);
        bool ok = true;
        for(int r=0; r<np; ++r)
        {
            if(all[3*r+2]!=all[2] || all[3*r+2]<1) ok = false;
            for(int i=all[3*r]; i<all[3*r]+all[3*r+1]; ++i)
            {
                if(i<0 || i>=gnx) { ok = false; break; }
                ++cover[i];
            }
        }
        for(int i=0; i<gnx && ok; ++i)
        if(cover[i]!=1) ok = false;
        if(!ok)
        {
            dist = false;
            return false;
        }

        const int NZ = all[2];
        comm = cm; dist = true; nxl = L.nx; i0 = gi0;
        cnt.assign(np,0); dsp.assign(np,0);
        std::vector<int> c6(np), d6(np);
        for(int r=0; r<np; ++r)
        {
            cnt[r] = all[3*r+1]*NZ; dsp[r] = all[3*r]*NZ;
            c6[r] = 6*cnt[r];       d6[r] = 6*dsp[r];
        }

        std::vector<double> loc, G((size_t)6*gnx*NZ);
        pack(L,prec,loc);
        MPI_Allgatherv(loc.data(),(int)loc.size(),MPI_DOUBLE,G.data(),c6.data(),d6.data(),MPI_DOUBLE,comm);
        factor_from(G,gnx,NZ);
        return true;
    }

    // u = A^-1 f on the interior of the level (f, u: the level arrays; several ranks: collective)
    void solve(const sc_level &L, const std::vector<double> &f, std::vector<double> &u) const
    {
        const long N = (long)nx*nz;
        std::vector<double> &y = work;
        y.resize((size_t)N);
        std::vector<double> &yl = workl;
        yl.resize((size_t)nxl*nz);
        for(int i=0; i<nxl; ++i)
        for(int k=0; k<nz; ++k)
        {
            const long q = L.idx(i,0,k);
            yl[(size_t)i*nz+k] = L.act[q] ? f[q] : 0.0;
        }
        if(dist)
        MPI_Allgatherv(yl.data(),(int)yl.size(),MPI_DOUBLE,y.data(),cnt.data(),dsp.data(),MPI_DOUBLE,comm);
        else
        y.swap(yl);

        for(long m=0; m<N; ++m)
        {
            const long c0 = std::max(0L,m-nz);
            const double *row = &a[(size_t)m*w + (c0-m+nz)];
            double s = y[m];
            for(long c=c0; c<m; ++c)
            s -= row[c-c0]*y[c];
            y[m] = s;
        }
        for(long m=N-1; m>=0; --m)
        {
            const long c1 = std::min(N-1,m+nz);
            const double *row = &a[(size_t)m*w + nz];
            double s = y[m];
            for(long c=m+1; c<=c1; ++c)
            s -= row[c-m]*y[c];
            y[m] = s/row[0];
        }

        const long o = (long)i0*nz;
        for(int i=0; i<nxl; ++i)
        for(int k=0; k<nz; ++k)
        u[L.idx(i,0,k)] = y[o + (long)i*nz+k];

        if(!dist)
        y.swap(yl);
    }

    mutable std::vector<double> work, workl;


private:

    // per unknown (i,k), k fastest: p, n, s, t, b, act of the level
    static void pack(const sc_level &L, int prec, std::vector<double> &G)
    {
        auto cf = [&](const std::vector<double> &d, const std::vector<float> &f, long q) -> double
        { return (prec==32) ? double(f[q]) : d[q]; };
        G.assign((size_t)6*L.nx*L.nz,0.0);
        for(int i=0; i<L.nx; ++i)
        for(int k=0; k<L.nz; ++k)
        {
            const long q = L.idx(i,0,k);
            double *g = &G[(size_t)6*((long)i*L.nz+k)];
            g[0] = cf(L.p,L.pf,q); g[1] = cf(L.n,L.nf,q); g[2] = cf(L.s,L.sf,q);
            g[3] = cf(L.t,L.tf,q); g[4] = cf(L.b,L.bf,q); g[5] = L.act[q] ? 1.0 : 0.0;
        }
    }

    void factor_from(const std::vector<double> &G, int NX, int NZ)
    {
        nx = NX; nz = NZ; w = 2*nz+1;
        const long N = (long)nx*nz;
        a.assign((size_t)N*w,0.0);
        auto act = [&](long m) { return G[(size_t)6*m+5]!=0.0; };

        for(int i=0; i<nx; ++i)
        for(int k=0; k<nz; ++k)
        {
            const long m = (long)i*nz + k;
            const double *g = &G[(size_t)6*m];
            double *row = &a[(size_t)m*w];
            if(!act(m))
            {
                row[nz] = 1.0;
                continue;
            }
            row[nz] = g[0];
            if(i+1<nx && act(m+nz)) row[nz+nz] = g[1];
            if(i-1>=0 && act(m-nz)) row[0]     = g[2];
            if(k+1<nz && act(m+1))  row[nz+1]  = g[3];
            if(k-1>=0 && act(m-1))  row[nz-1]  = g[4];
        }

        // elimination: row r, column c at a[r*w + (c-r+nz)]
        for(long m=0; m<N; ++m)
        {
            const double piv = a[(size_t)m*w+nz];
            const long r1 = std::min(N-1,m+nz);
            for(long r=m+1; r<=r1; ++r)
            {
                double &lrm = a[(size_t)r*w + (m-r+nz)];
                if(lrm==0.0)
                continue;
                lrm /= piv;
                const long c1 = std::min(N-1,m+nz);
                const double *um = &a[(size_t)m*w + nz];      // um[c-m]: U(m,c)
                double *rr = &a[(size_t)r*w + (m-r+nz)];      // rr[c-m]: A(r,c)
                for(long c=m+1; c<=c1; ++c)
                rr[c-m] -= lrm*um[c-m];
            }
        }

    }
};

#endif
