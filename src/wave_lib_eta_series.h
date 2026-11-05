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

#ifndef WAVE_LIB_ETA_SERIES_H_
#define WAVE_LIB_ETA_SERIES_H_

#include<vector>
#include<cmath>

/*--------------------------------------------------------------------
Time series of a sum of harmonics at one point (iowave::timeseries,
REEF3D_Log-Wave):

  eta(t) = sum_b ( ar_b cos(W_b t) + ai_b sin(W_b t) )

From one sample to the next the harmonics are rotated by W_b dt (complex
multiplication, no cos / sin per sample); they are evaluated exactly every
32 samples and whenever the step differs from the first one by more than
1e-9 of it. Agrees with the direct evaluation to about 1e-12 of the
amplitude sum.
--------------------------------------------------------------------*/

inline void wave_eta_harmonics(const std::vector<double> &W, const std::vector<double> &ar, const std::vector<double> &ai,
                               const std::vector<double> &tv, std::vector<double> &ev)
{
    const int B = int(W.size());
    const size_t NT = tv.size();
    ev.assign(NT,0.0);

    if(B==0 || NT==0)
    return;

    std::vector<double> c(B), s(B), cd(B,1.0), sd(B,0.0);
    const double d0 = NT>1 ? tv[1] - tv[0] : 0.0;

    for(int b=0; b<B; ++b)
    {
    cd[b] = cos(W[b]*d0);
    sd[b] = sin(W[b]*d0);
    }

    const double *pw = W.data(), *pr = ar.data(), *pi = ai.data(), *pcd = cd.data(), *psd = sd.data();
    double *pc = c.data(), *ps = s.data();
    int since = 32;

    for(size_t q=0; q<NT; ++q)
    {
        const double t = tv[q];
        const bool exact = since>=32 || q==0 || fabs((tv[q] - tv[q-1]) - d0) > 1.0e-9*fabs(d0);

        if(exact)
        {
            for(int b=0; b<B; ++b)
            {
            pc[b] = cos(pw[b]*t);
            ps[b] = sin(pw[b]*t);
            }
            since = 1;
        }
        else
        {
            #pragma omp simd
            for(int b=0; b<B; ++b)
            {
                const double cn = pc[b]*pcd[b] - ps[b]*psd[b];
                ps[b] = ps[b]*pcd[b] + pc[b]*psd[b];
                pc[b] = cn;
            }
            ++since;
        }

        double e = 0.0;

        #pragma omp simd reduction(+:e)
        for(int b=0; b<B; ++b)
        e += pr[b]*pc[b] + pi[b]*ps[b];

        ev[q] = e;
    }
}

#endif
