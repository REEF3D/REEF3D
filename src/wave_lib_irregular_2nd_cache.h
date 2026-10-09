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

#ifndef WAVE_LIB_IRREGULAR_2ND_CACHE_H_
#define WAVE_LIB_IRREGULAR_2ND_CACHE_H_

#include<vector>
#include<cmath>

/*--------------------------------------------------------------------
Cached-point evaluation of the 2nd-order irregular theories
(wave_lib_irregular_2nd_a / _b, B 92 32, 33, 42, 43, 52, 53).

The direct evaluation costs cosh / sinh / cos of every component pair at
every point and call, O(M^2) transcendental functions for M components.
Here the pair terms are built from per-component quantities:

  phases     C_n = cos T_n, S_n = sin T_n with T_n = a_n(x,y) + b_n(t):
             cos a_n, sin a_n per registered point (once), cos b_n,
             sin b_n per time (once per step), angle addition per call;
             cos(T_n -+ T_m) = C_n C_m +- S_n S_m, sin(T_n -+ T_m) = ...
  vertical   E_n = exp(k_n s), R_n = 1 / E_n with s = d + z, per call;
             cosh((k_n -+ k_m) s) = (E_n R_m + R_n E_m) / 2 or
             (E_n E_m + R_n R_m) / 2, sinh likewise

so a call costs M exp and O(M^2) multiplications. Same terms as the
direct evaluation, results equal to round-off. When k_max s > 340 the
products could overflow and the library falls back to the direct
evaluation for that call (very deep water only).
--------------------------------------------------------------------*/

struct wave_lib_irregular_2nd_cache
{
    int M=0;                              // components
    std::vector<double> cS, sS;           // cos / sin of the spatial phase, point-major
    std::vector<double> cT, sT;           // cos / sin of the temporal phase
    double t=-1.0e300;
    std::vector<double> C, S, E, R;       // per call
    std::vector<double> k;                // wave numbers
    double kmax=0.0;
    bool on=false;

    void points(const std::vector<double> &x, const std::vector<double> &y, int m,
                const double *ki, const double *cosb, const double *sinb)
    {
        M = m;
        const int N = int(x.size());
        cS.assign(size_t(N)*M,0.0);
        sS.assign(size_t(N)*M,0.0);

        for(int q=0; q<N; ++q)
        for(int n=0; n<M; ++n)
        {
            const double a = ki[n]*(cosb[n]*x[q] + sinb[n]*y[q]);
            cS[size_t(q)*M+n] = cos(a);
            sS[size_t(q)*M+n] = sin(a);
        }

        cT.assign(M,0.0); sT.assign(M,0.0);
        C.assign(M,0.0); S.assign(M,0.0); E.assign(M,0.0); R.assign(M,0.0);
        k.assign(ki,ki+M);
        kmax = 0.0;
        for(int n=0; n<M; ++n)
        kmax = fmax(kmax,fabs(ki[n]));

        t = -1.0e300;
        on = true;
    }

    void time(double wt, const double *wi, const double *ei)
    {
        if(wt==t)
        return;

        for(int n=0; n<M; ++n)
        {
            const double b = -wi[n]*wt - ei[n];
            cT[n] = cos(b);
            sT[n] = sin(b);
        }
        t = wt;
    }

    // phases of point q at the cached time
    void phases(int q)
    {
        const double *cs=&cS[size_t(q)*M], *ss=&sS[size_t(q)*M];

        for(int n=0; n<M; ++n)
        {
            C[n] = cs[n]*cT[n] - ss[n]*sT[n];
            S[n] = ss[n]*cT[n] + cs[n]*sT[n];
        }
    }

    // vertical factors at s = d + z; false: use the direct evaluation
    bool vertical(double s)
    {
        if(kmax*fabs(s)>340.0)
        return false;

        for(int n=0; n<M; ++n)
        {
            E[n] = exp(k[n]*s);
            R[n] = 1.0/E[n];
        }
        return true;
    }
};

#endif
