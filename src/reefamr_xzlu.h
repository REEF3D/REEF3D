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

//  Exact solve of the 7-point operator of a reefmg_core fine level with one cell in y (a patch of a
//  2D grid: an x-z problem), in place of its V-cycle.  The unknowns are ordered column by column
//  (k fastest), so the matrix is banded with nz sub- and super-diagonals: LU without pivoting in
//  band storage (the rows are those the line smoother of the V-cycle factorises without pivoting
//  too), O(nx nz^2) per factorisation, O(nx nz) per solve.  The couplings beyond the level (the
//  patch edge, below the bed node, above the top node) are 0 in the rows or left out, rows with
//  act 0 are identity rows.
struct reefamr_xzlu
{
    int nx = 0, nz = 0, w = 0;         // w = 2 nz + 1, the band width
    std::vector<double> a;              // row m: entries of the columns m-nz .. m+nz

    bool usable(const sc_level &L) const { return L.ny==1; }

    void factor(const sc_level &L, int prec)
    {
        nx = L.nx; nz = L.nz; w = 2*nz+1;
        const long N = (long)nx*nz;
        a.assign((size_t)N*w,0.0);

        auto cf = [&](const std::vector<double> &d, const std::vector<float> &f, long q) -> double
        { return (prec==32) ? double(f[q]) : d[q]; };

        for(int i=0; i<nx; ++i)
        for(int k=0; k<nz; ++k)
        {
            const long m = (long)i*nz + k;
            const long q = L.idx(i,0,k);
            double *row = &a[(size_t)m*w];
            if(L.act[q]==0)
            {
                row[nz] = 1.0;
                continue;
            }
            row[nz] = cf(L.p,L.pf,q);
            if(i+1<nx && L.act[L.idx(i+1,0,k)]) row[nz+nz] = cf(L.n,L.nf,q);
            if(i-1>=0 && L.act[L.idx(i-1,0,k)]) row[0]     = cf(L.s,L.sf,q);
            if(k+1<nz && L.act[q+1])            row[nz+1]  = cf(L.t,L.tf,q);
            if(k-1>=0 && L.act[q-1])            row[nz-1]  = cf(L.b,L.bf,q);
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

    // u = A^-1 f on the interior of the level (f, u: the level arrays)
    void solve(const sc_level &L, const std::vector<double> &f, std::vector<double> &u) const
    {
        const long N = (long)nx*nz;
        std::vector<double> &y = work;
        y.resize((size_t)N);
        for(int i=0; i<nx; ++i)
        for(int k=0; k<nz; ++k)
        {
            const long q = L.idx(i,0,k);
            y[(size_t)i*nz+k] = L.act[q] ? f[q] : 0.0;
        }

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

        for(int i=0; i<nx; ++i)
        for(int k=0; k<nz; ++k)
        u[L.idx(i,0,k)] = y[(size_t)i*nz+k];
    }

    mutable std::vector<double> work;
};

#endif
