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

#ifndef FNPF_SZ_WEIGHTS_H_
#define FNPF_SZ_WEIGHTS_H_

#include"lexer.h"
#include"increment.h"

/*
    One-sided d/dsigma at the free-surface node k = knoz (A327 1).

    The weights are computed on the actual sigma-node positions ZN
    (Fornberg 1988, Math. Comp. 51, 699-706), so an n-point stencil is exact
    for polynomials of degree n-1 in sigma for any node spacing.

    The legacy form sum(c*f)/sum(c*ZN) with uniform-grid coefficients c
    (A327 0) is the chain rule in grid-index space: exact for polynomials in
    the index, not in sigma. It is accurate only while ZN(k) is smooth on the
    scale of one index step and degrades on strongly stretched vertical grids,
    e.g. few layers in deep water (B103 5, large B113).

    ZN is the same for all columns, so the weights are computed once; the
    caller converts to d/dz with sigz.

    Stencils: 7 (cds6), 5 (cds4, weno5), 3 (cds2). If a node of the stencil
    is not active (flag7<=0) or there are too few layers, the next smaller
    stencil 5 -> 3 -> 2 is used.
*/

class fnpf_sz_weights
{
public:
    fnpf_sz_weights(lexer *p)
    {
        const int kt = p->knoz;
        const int mg = increment::marge;

        for(int n=0; n<=maxpt; ++n)
        {
            valid[n] = false;
            for(int m=0; m<maxpt; ++m)
            w[n][m] = 0.0;
        }

        for(int n=2; n<=maxpt; ++n)
        if(kt>=n-1)
        {
            double x[maxpt];

            for(int m=0; m<n; ++m)
            x[m] = p->ZN[kt-m+mg];   // x[0]: surface node, x[m]: node k-m

            ddz_weights(x[0],x,n,w[n]);
            valid[n] = true;
        }
    }

    // f: FNPF potential (c->Fi), q: FIJK of the surface node, npt: stencil points
    inline double eval(lexer *p, const double *f, int q, int npt) const
    {
        static const int seq[4] = {7,5,3,2};

        for(int s=0; s<4; ++s)
        {
            const int n = seq[s];

            if(n>npt || !valid[n])
            continue;

            bool active = true;
            for(int m=0; m<n-1; ++m)   // as legacy: lowest stencil node not checked
            if(p->flag7[q-m]<=0)
            {
                active = false;
                break;
            }

            if(active || n==2)
            {
                double val = 0.0;
                for(int m=0; m<n; ++m)
                val += w[n][m]*f[q-m];

                return val;
            }
        }

        return 0.0;
    }

private:
    static const int maxpt = 7;

    // Fornberg (1988): weights of the first derivative at z0 on the nodes x[0..n-1]
    static void ddz_weights(double z0, const double *x, int n, double *out)
    {
        double c[maxpt][2];

        for(int r=0; r<n; ++r)
        c[r][0] = c[r][1] = 0.0;

        double c1 = 1.0;
        double c4 = x[0] - z0;
        c[0][0] = 1.0;

        for(int r=1; r<n; ++r)
        {
            const int mn = r<1?r:1;
            double c2 = 1.0;
            const double c5 = c4;
            c4 = x[r] - z0;

            for(int s=0; s<r; ++s)
            {
                const double c3 = x[r] - x[s];
                c2 *= c3;

                if(s==r-1)
                {
                    for(int d=mn; d>=1; --d)
                    c[r][d] = c1*(double(d)*c[r-1][d-1] - c5*c[r-1][d])/c2;

                    c[r][0] = -c1*c5*c[r-1][0]/c2;
                }

                for(int d=mn; d>=1; --d)
                c[s][d] = (c4*c[s][d] - double(d)*c[s][d-1])/c3;

                c[s][0] = c4*c[s][0]/c3;
            }
            c1 = c2;
        }

        for(int r=0; r<n; ++r)
        out[r] = c[r][1];
    }

    double w[maxpt+1][maxpt];
    bool valid[maxpt+1];
};

#endif
