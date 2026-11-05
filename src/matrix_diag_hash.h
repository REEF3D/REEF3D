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

#ifndef MATRIX_DIAG_HASH_H_
#define MATRIX_DIAG_HASH_H_

#include"matrix_diag.h"
#include<cstring>
#include<vector>

// 64-bit hash (FNV-1a on the bit patterns, four lanes) of the 7-point rows of a matrix_diag, the
// cross-couplings nt..eb included when they are sized: a solver whose coarse levels or
// factorisations are made from the rows keeps them for the next solve with the same hash (REEFMG:
// the psi solves of the FNPF floating bodies have the rows of the phi solve before them).  Only
// what is made from the rows for the preconditioner may be kept this way: a collision (2^-64) would
// cost convergence, not accuracy, as long as the Krylov operator reads the rows themselves
struct matrix_diag_hash_state
{
    unsigned long long h[4] = {1469598103934665603ULL, 0x84222325CBF29CE4ULL, 0x9E3779B97F4A7C15ULL, 0xC2B2AE3D27D4EB4FULL};
    
    void add(const std::vector<double> &v)
    {
        const unsigned long long P = 1099511628211ULL;
        const size_t n = v.size();
        size_t m = 0;
        unsigned long long w[4];
        
        for(; m+4<=n; m+=4)
        {
            std::memcpy(w,&v[m],4*sizeof(double));
            h[0] = (h[0]^w[0])*P;
            h[1] = (h[1]^w[1])*P;
            h[2] = (h[2]^w[2])*P;
            h[3] = (h[3]^w[3])*P;
        }
        
        for(; m<n; ++m)
        {
            std::memcpy(w,&v[m],sizeof(double));
            h[0] = (h[0]^w[0])*P;
        }
        
        h[1] = (h[1]^(unsigned long long)n)*P;
    }
    
    void add(const matrix_diag &M)
    {
        add(M.p); add(M.n); add(M.s); add(M.w); add(M.e); add(M.t); add(M.b);
        add(M.nt); add(M.nb); add(M.st); add(M.sb); add(M.wt); add(M.wb); add(M.et); add(M.eb);
    }
    
    unsigned long long value() const
    {
        return h[0] ^ (h[1]<<1 | h[1]>>63) ^ (h[2]<<2 | h[2]>>62) ^ (h[3]<<3 | h[3]>>61);
    }
};

#endif
