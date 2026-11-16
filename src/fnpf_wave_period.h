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

#ifndef FNPF_WAVE_PERIOD_H_
#define FNPF_WAVE_PERIOD_H_

#include"lexer.h"
#include<fstream>
#include<sstream>
#include<string>
#include<vector>
#include<algorithm>
#include<cmath>

// The period of the wave input for the relaxation zones of theta-Robin (A 303 <= -1, fnpf_theta::relax_scale):
// regular waves wT, spectra wTp (as the step of N 48 1, fnpf_timestep::theta_tref, which does not use the peak
// period of reconstructed waves). Reconstructed waves (B 92 51-53,
// waverecon.dat "A omega phase [direction]") have no wTp (other code takes wTp as the input of a spectrum, so it
// stays 0 there); their period is the peak of the component spectrum S_n = A_n^2/(2 dw_n), smoothed as
// Fp = sum(w S^4)/sum(S^4) (Young 1995), dw_n half the distance of the neighbouring frequencies. 0 without wave
// input.

inline double fnpf_recon_period()
{
    std::ifstream file("waverecon.dat", std::ios_base::in);

    if(!file)
    return 0.0;

    std::vector<double> A, w;
    std::string line;

    while(std::getline(file,line))
    {
        std::istringstream ls(line);
        double a, om;

        if(ls>>a>>om)
        if(om>0.0)
        {
        A.push_back(a);
        w.push_back(om);
        }
    }

    const int N = int(w.size());

    if(N<1)
    return 0.0;

    std::vector<int> o(N);
    for(int n=0;n<N;++n)
    o[n]=n;
    std::sort(o.begin(),o.end(),[&](int a, int b){return w[a]<w[b];});

    double s4=0.0, ws4=0.0;
    for(int m=0;m<N;++m)
    {
        const int n = o[m];
        double dw = 1.0;

        if(N>1)
        dw = (m==0) ? w[o[1]]-w[o[0]] : ((m==N-1) ? w[o[N-1]]-w[o[N-2]] : 0.5*(w[o[m+1]]-w[o[m-1]]));

        if(dw<=0.0)
        continue;

        const double S = A[n]*A[n]/(2.0*dw);
        const double S4 = S*S*S*S;
        s4 += S4;
        ws4 += w[n]*S4;
    }

    return (s4>0.0 && ws4>0.0) ? 2.0*M_PI*s4/ws4 : 0.0;
}

inline double fnpf_wave_period(lexer *p)
{
    if(p->B90<=0)
    return 0.0;

    if(p->B92<=11)
    return p->wT;

    if(p->wTp>0.0)
    return p->wTp;

    if(p->B92>=51 && p->B92<=53)
    {
        static double T = -1.0;

        if(T<0.0)
        T = fnpf_recon_period();

        return T;
    }

    return 0.0;
}

#endif
