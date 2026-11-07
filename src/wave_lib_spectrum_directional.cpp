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
Author: Widar W. Wang
--------------------------------------------------------------------*/

#include"wave_lib_spectrum.h"
#include"lexer.h"
#include"ghostcell.h"
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand
#include <random>
#include <algorithm>
#include <numeric>
#include <vector>

void wave_lib_spectrum::directional_spreading(lexer* p) // modified
{
    int n, q;
    int count;

    // calc delta beta
    if(p->B130==0)
    dbeta = 1.0;

    if(p->B130>0)
    {

        /*if(p->B133 % 2==0)
            p->B133 += 1;*/

        if(p->B133 <= 1)
            p->B133 = 3;

        // convert beta_start and beta_end and calculate delta_beta
        p->B131 *= PI / 180.0;
        p->B132_s *= PI / 180.0;
        p->B132_e *= PI / 180.0;

        // =================================================================
        // Special handling for 2D spectrum file (B85==11)
        // Use file frequencies and directions directly - no interpolation
        // =================================================================
        if(p->B85==11)
        {
            // Use file directions directly (p->B133 already set in wave_lib_spectrum.cpp)

            // Setup beta array using file directions
            count = 0;
            for(n = 0; n < p->wN; ++n)
            {
                for(q = 0; q < p->B133; ++q)
                {
                    beta[count] = dir_2d[q];
                    sinbeta[count] = sin(beta[count]);
                    cosbeta[count] = cos(beta[count]);
                    ++count;
                }
            }

            // Calculate dbeta (average directional spacing from file)
            if(ptnum_dir_2d > 1)
                dbeta = (dir_2d[ptnum_dir_2d-1] - dir_2d[0]) / double(ptnum_dir_2d - 1);
            else
                dbeta = 0.0;

            // Setup wi array - expand frequencies
            double* Si_temp;
            p->Darray(Si_temp, p->wN);

            for(n = 0; n < p->wN; ++n)
                Si_temp[n] = wi[n];

            count = 0;
            for(n = 0; n < p->wN; ++n)
            {
                for(q = 0; q < p->B133; ++q)
                {
                    wi[count] = Si_temp[n];
                    ++count;
                }
            }

            // Setup dw array - expand frequency intervals
            for(n = 0; n < p->wN; ++n)
                Si_temp[n] = dw[n];

            count = 0;
            for(n = 0; n < p->wN; ++n)
            {
                for(q = 0; q < p->B133; ++q)
                {
                    dw[count] = Si_temp[n];
                    ++count;
                }
            }

            // Setup ki array - expand to freq*dir components
            for(n = 0; n < p->wN; ++n)
                Si_temp[n] = ki[n];

            count = 0;
            for(n = 0; n < p->wN; ++n)
            {
                for(q = 0; q < p->B133; ++q)
                {
                    ki[count] = Si_temp[n];
                    ++count;
                }
            }

            // Get spectrum values directly from file data (no interpolation)
            count = 0;
            for(n = 0; n < p->wN; ++n)
            {
                for(q = 0; q < p->B133; ++q)
                {
                    Di[count] = 1.0;

                    // Direct indexing into spectrum_2d_1d array
                    Si[count] = spectrum_2d_1d[n * ptnum_dir_2d + q];

                    ++count;
                }
            }

            p->del_Darray(Si_temp, p->wN);

            // Setup beta_n and Di_n arrays for print_spreading
            for(q = 0; q < p->B133; ++q)
            {
                beta_n[q] = dir_2d[q];
                Di_n[q] = 1.0;  // Flat spreading since file already contains directional info
            }

            // Update p->wN to total components
            p->wN *= p->B133;

            if(p->mpirank==0)
                cout<<"2D Spectrum: dbeta = "<<dbeta<<" rad"<<endl;
        }
        // =================================================================
        // End of B85==11 special handling
        // Below: Standard directional spreading (analytical methods)
        // =================================================================
        else if(p->B136==4 && p->B84!=2)
        {
            if(p->mpirank==0)
            cout<<"!!! B 136 4 (equal energy directions) needs B 84 2 (equal energy frequencies) !!!"<<endl;
            exit(1);
        }

        // B 136 1, 2, 3: B 133 directions of equal width per frequency, also for B 84 2 (it gave unidirectional
        // waves at 0 deg). The sector B 132 is relative to B 131 here
        else if(p->B136!=4)
        {
            dbeta = (p->B132_e - p->B132_s) / double(p->B133 - 1);

            if(p->mpirank==0)
            cout<<"directional spreading: B 136 "<<p->B136<<", sector "<<(p->B131+p->B132_s)*180.0/PI<<" to "<<(p->B131+p->B132_e)*180.0/PI
                <<" deg (B 132 relative to B 131 "<<p->B131*180.0/PI<<" deg)"<<endl;

            // calc beta[i]
            count = 0;
            for(n = 0; n < p->wN; ++n)
            {
                for(q = 0; q < p->B133; ++q)
                {
                    beta[count] = p->B132_s + dbeta * double(q) + p->B131;
                    sinbeta[count] = sin(beta[count]);
                    cosbeta[count] = cos(beta[count]);
                    ++count;
                }
            }

            double* Si_temp;
            p->Darray(Si_temp, p->wN);

            // Calc wi
            for(n = 0; n < p->wN; ++n)
                Si_temp[n] = wi[n];

            count = 0;
            // Double Summation Method frequency vector
            if (p->B136==1)
            {
                for(n = 0; n < p->wN; ++n)
                {
                    for(q = 0; q < p->B133; ++q)
                    {
                        wi[count] = Si_temp[n];
                        ++count;
                    }
                }
            }
            // Pascal's Single Summation Method frequency vector
            if (p->B136==2)
            {
                for(n = 0; n < p->wN; ++n)
                {
                    for(q = 0; q < p->B133; ++q)
                    {
                        wi[count] = Si_temp[n]+double (q)*dw[n]/double (p->B133);
                        ++count;
                    }
                }
            }

            // Yu's Single Summation Method frequency vector
            if (p->B136==3)
            {
                if(p->B138>0)
                srand(p->B138_1);

                if(p->B138==0)
                srand((unsigned)time(0));
                
                for(n = 0; n < p->wN; ++n)
                {
                    for(q = 0; q < p->B133; ++q)
                    {
                        wi[count] = Si_temp[n]+(double (q)+((float) rand()) / (float) RAND_MAX)*dw[n]/double (p->B133);
                        ++count;
                    }
                }
            }

            // Calc dw
            for(n = 0; n < p->wN; ++n)
            Si_temp[n] = dw[n];

            count = 0;
            for(n = 0; n < p->wN; ++n)
            {
                for(q = 0; q < p->B133; ++q)
                {
                dw[count] = Si_temp[n];
                ++count;
                }
            }

            // Calc ki
            for(n = 0; n < p->wN; ++n)
                Si_temp[n] = ki[n];

            count = 0;
            for(n = 0; n < p->wN; ++n)
            {
                for(q = 0; q < p->B133; ++q)
                {
                    ki[count] = Si_temp[n];
                    ++count;
                }
            }

            // single summation (B 136 2, 3): every component has its own frequency, so its own
            // wave number; k of the frequency bin gave the wrong wave length in the generation zone
            if(p->B136==2 || p->B136==3)
            for(n = 0; n < p->wN*p->B133; ++n)
            {
                const double w2 = wi[n]*wi[n];
                double kk = ki[n];

                for(int it = 0; it < 20; ++it)
                {
                    const double kd = kk*wD;
                    const double dk = (9.81*kk*tanh(kd) - w2)/(9.81*(tanh(kd) + kd/(cosh(kd)*cosh(kd))));
                    kk -= dk;

                    if(fabs(dk) < 1.0e-12*kk)
                    break;
                }

                ki[n] = kk;
            }

            for(n = 0; n < p->wN; ++n)
                Si_temp[n] = Si[n];

            // multiply S with D(beta), using spreading function. D is scaled per frequency to
            // sum(D dbeta) = 1 over the discrete directions: the energy of the frequency stays in
            // the sector, whatever the shape parameter, B 133 and B 132 (the form of B 130 1 has
            // the integral 1 only for s = 2, with B 134 10 the waves had 0.70 of the Hs)
            count = 0;
            for(n = 0; n < p->wN; ++n)
            {
                double sumD = 0.0;

                for(q = 0; q < p->B133; ++q)
                {
                  Di[count+q] = spreading_function(p, beta[count+q], wi[count+q]);
                  sumD += Di[count+q]*dbeta;
                }

                if(n==0 && p->mpirank==0 && sumD < 0.9)
                cout<<"!!! the sector B 132 holds "<<100.0*sumD<<" % of the spreading function at the first frequency !!!"<<endl;

                for(q = 0; q < p->B133; ++q)
                {
                  if(sumD>0.0)
                  Di[count] /= sumD;

                  Si[count] = Si_temp[n] * Di[count];
                  ++count;
                }
            }

            p->del_Darray(Si_temp, p->wN);

            // set p->wN to B86*B133
            p->wN *= p->B133;

            if(p->mpirank==0)
            {
                cout << "directional spreading  dbeta: " << dbeta << "  Nj: " << p->B133 << endl;
            }
        }

        else if(p->B84==2 && p->B136==4)
        {
            // equal energy directions: every component has the energy of its frequency (equal energy, B 84 2) and
            // takes the direction at the middle of one of B 133 intervals of equal probability of D(beta|w), the
            // distribution of its own frequency, so a spreading that changes with frequency (B 135) is kept.
            // Each block of B 133 neighbouring frequencies holds every interval once, in a random order. The
            // sector B 132 is absolute here (in the frame of the model), as before.
            // Before: the directions of the levels 0.5 % ... 99.5 % of the cdf, the end levels included (spreading
            // about 6 % too wide), one D(beta) for all frequencies with the frequency running along the
            // directions (B 135 gave the wrong spreading and mean), and a shuffle over all frequencies.
            const int Nd = p->B133, N = p->wN;

            if(Nd<1 || Nd>N || N % Nd!=0)
            {
                if(p->mpirank==0)
                cout<<"!!! B 133 (number of directions "<<Nd<<") must divide the number of components "<<N<<": the next ones are "
                    <<(N/std::max(Nd,1))*Nd<<" and "<<(N/std::max(Nd,1)+1)*Nd<<" for B 86 !!!"<<endl;
                exit(1);
            }

            const double ddb = 0.002;
            const int ND = std::max(int((p->B132_e - p->B132_s) / ddb) + 1, 2);
            const double hb = (p->B132_e - p->B132_s)/double(ND-1);
            const bool perfreq = (p->B130==2 && p->B134<=0.0);

            std::vector<double> bt(ND), cdfb(ND,0.0);

            for(int m=0; m<ND; ++m)
            bt[m] = p->B132_s + double(m)*hb;

            auto build = [&](const double w)
            {
                double prev = spreading_function(p,bt[0],w);
                cdfb[0] = 0.0;

                for(int m=1; m<ND; ++m)
                {
                    const double v = spreading_function(p,bt[m],w);
                    cdfb[m] = cdfb[m-1] + 0.5*(prev+v)*hb;
                    prev = v;
                }
            };

            auto invert = [&](const double u)
            {
                const double v = u*cdfb[ND-1];
                int m = int(std::lower_bound(cdfb.begin(),cdfb.end(),v) - cdfb.begin());
                m = std::min(std::max(m,1),ND-1);
                const double h = cdfb[m]-cdfb[m-1];

                return h>0.0 ? bt[m-1] + (bt[m]-bt[m-1])*(v-cdfb[m-1])/h : bt[m];
            };

            // the levels at the peak frequency, for the log and REEF3D_spreading-function.dat
            build(p->wwp);

            if(p->mpirank==0)
            {
                cout<<"directional spreading: B 136 4, "<<Nd<<" equal energy directions, sector "<<p->B132_s*180.0/PI<<" to "<<p->B132_e*180.0/PI
                    <<" deg (B 132 absolute), mean B 131 "<<p->B131*180.0/PI<<" deg"<<endl;

                if(cdfb[ND-1] < 0.9)
                cout<<"!!! the sector B 132 holds "<<100.0*cdfb[ND-1]<<" % of the spreading function: with B 136 4, B 132 is in the frame of the model, not relative to B 131 !!!"<<endl;
            }

            for(int q=0; q<Nd; ++q)
            {
                beta_n[q] = invert((double(q)+0.5)/double(Nd));
                Di_n[q] = spreading_function(p,beta_n[q],p->wwp);
            }

            std::mt19937_64 gen(p->B138>0 ? (unsigned long long)p->B138_2 : (unsigned long long)time(0));
            std::vector<int> lvl(Nd);
            std::iota(lvl.begin(),lvl.end(),0);

            if(!perfreq)
            build(p->wwp);

            for(n = 0; n < N; ++n)
            {
                if(n % Nd==0)
                std::shuffle(lvl.begin(),lvl.end(),gen);

                if(perfreq)
                build(wi[n]);

                beta[n] = invert((double(lvl[n % Nd])+0.5)/double(Nd));
                sinbeta[n] = sin(beta[n]);
                cosbeta[n] = cos(beta[n]);

                // Si for the component file; the amplitudes use Sn (equal energy)
                Di[n] = spreading_function(p,beta[n],wi[n]);
                Si[n] = Sn[n]*Di[n];
            }
        }
    }
}

double wave_lib_spectrum::spreading_function(lexer* p, double beta, double w)
{
    // D(beta) normalised to 1 over the circle, the difference to the mean direction B 131 taken in (-pi, pi],
    // so that a sector given past -180 or 180 deg still finds the peak
    double D = 0.0;
    const double d = atan2(sin(beta - p->B131), cos(beta - p->B131));

    // B 130 1: cos^s(beta - B 131) within 90 deg of B 131, the form of SWAN (DSPR POWER, m = B 134);
    // any s > 0. The factor (2/pi) was the normalisation of s = 2 only
    if(p->B130==1)
    {
        s_f = p->B134;

        if(fabs(d) < 0.5*PI && s_f > 0.0)
        D = exp(lgamma(0.5*s_f + 1.0) - lgamma(0.5*s_f + 0.5))/sqrt(PI) * pow(cos(d), s_f);
    }

    // B 130 2: cos^2s((beta - B 131)/2) over the circle (Mitsuyasu, Longuet-Higgins). B 134 > 0 gives s, otherwise
    // s of the frequency from the peak value B 135 (Goda: (w/wp)^5 below the peak, (w/wp)^-2.5 above it).
    // G0 from log-gamma: Gamma(2s+1) overflowed above s = 85, which gave no waves (NaN)
    if(p->B130==2)
    {
        if(p->B134 > 0.0)
        s_f = p->B134;

        else
        s_f = (w <= p->wwp) ? p->B135 * pow(w / p->wwp, 5.0) : p->B135 * pow(w / p->wwp, -2.5);

        G_0 = exp((2.0*s_f - 1.0)*log(2.0) - log(PI) + 2.0*lgamma(s_f + 1.0) - lgamma(2.0*s_f + 1.0));

        D = G_0 * pow(fabs(cos(0.5*d)), 2.0 * s_f);
    }

    return D;
}
