/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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
Author: Weizhi Wang
--------------------------------------------------------------------*/

#include "wave_lib_spectrum.h"
#include "lexer.h"
#include "ghostcell.h"
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand

void wave_lib_spectrum::directional_spreading(lexer* p) // modified
{
    int n, q;
    int count;

    // calc delta beta
    if(p->B130 == 0)
    dbeta = 1.0;

    if(p->B130>0)
    {

        /*if(p->B133 % 2 == 0)
            p->B133 += 1;*/

        if(p->B133 <= 1)
            p->B133 = 3;

        // convert beta_start and beta_end and calculate delta_beta
        p->B131 *= PI / 180.0;
        p->B132_s *= PI / 180.0;
        p->B132_e *= PI / 180.0;

        if(p->B84==1)
        {
            dbeta = (p->B132_e - p->B132_s) / double(p->B133 - 1);

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
            if (p->B136 == 1)
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
            if (p->B136 == 2)
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
            if (p->B136 == 3)
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

            for(n = 0; n < p->wN; ++n)
                Si_temp[n] = Si[n];

            // multiply S with D(beta), depending on spreading function
            count = 0;
            for(n = 0; n < p->wN; ++n)
            {
                for(q = 0; q < p->B133; ++q)
                {
                  Di[count]=spreading_function(p, beta[count], wi[count]);
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

        if(p->B84==2 && p->B136==4)
        {
            double DD, ddb, w, ddw, sum, Dcdf_s, Dcdf_e;
            double cdf_s, cdf_e, d_s, d_e, d_low, d_high, cdf_low, cdf_high;
            int m, NR, ND;

            ddb=0.01;
            ND= (int) ((p->B132_e - p->B132_s) / ddb)+1;

            // integration of the spectrum
            p->Darray(Dn,ND);
            p->Darray(Dcdf,ND);
            p->Darray(betat,ND);
            p->Darray(Ddee,p->B133);

            w=ws;
            ddw=(we-ws)/double(ND-1);
            for(q = 0; q < ND; ++q)
            {
                betat[q]=p->B132_s+double (q)*ddb;
                DD=spreading_function(p,betat[q],w);
                Dn[q]=DD;
                w+=ddw;
            }

            sum=0.0;
            for(n=0;n<ND;++n)
            {
                sum+=ddb*Dn[n];
                Dcdf[n]=sum;
            }

            // Create 0.5% - 99.5% caps of the energy
            Dcdf_s = 0.005*Dcdf[ND-1];
            Dcdf_e = 0.995*Dcdf[ND-1];

            for(m=0;m<ND;++m)
            {
                if(Dcdf[m]<=Dcdf_s)
                {
                    cdf_low=Dcdf[m];
                    d_low=betat[m];
                }
            }

            for(m=(ND-1);m>=0;--m)
            {
                if(Dcdf[m]>=Dcdf_s)
                {
                    cdf_high=Dcdf[m];
                    d_high=betat[m];
                }
            }

            if(d_low==d_high)
            {
                d_s=d_low;
            }

            if(d_low!=d_high)
            {
                d_s=(Dcdf_s-cdf_low)*(d_high-d_low)/(cdf_high-cdf_low)+d_low;
            }

            for(m=0;m<ND;++m)
            {
                if(Dcdf[m]<=Dcdf_e)
                {
                    cdf_low=Dcdf[m];
                    d_low=betat[m];
                }
            }

            for(m=(ND-1);m>=0;--m)
            {
                if(Dcdf[m]>=Dcdf_e)
                {
                    cdf_high=Dcdf[m];
                    d_high=betat[m];
                }
            }

            if(d_low==d_high)
            {
                d_e=d_high;
            }

            if(d_low!=d_high)
            {
                d_e=(Dcdf_e-cdf_low)*(d_high-d_low)/(cdf_high-cdf_low)+d_low;
            }

            // Create equal energy bins cdf_s:(cdf_e-cdf_s)/double (p->B133-1):cdf_e
            for(n=0;n<p->B133;++n)
            {
                Ddee[n] = Dcdf_s+double (n)*(Dcdf_e-Dcdf_s)/double (p->B133-1);
            }

            // Interpolate the corresponding angles at each equal energy bin
            for(n=0;n<p->B133;++n)
            {
                for(m=0;m<ND;++m)
                {
                    if(Dcdf[m]<=Ddee[n])
                    {
                        cdf_low=Dcdf[m];
                        d_low=betat[m];
                    }
                }

                for(m=(ND-1);m>=0;--m)
                {
                    if(Dcdf[m]>=Ddee[n])
                    {
                        cdf_high=Dcdf[m];
                        d_high=betat[m];
                    }
                }

                if(d_low==d_high)
                {
                    beta[n]=d_low;
                }

                if(d_low!=d_high)
                {
                    beta[n]=(Ddee[n]-cdf_low)*(d_high-d_low)/(cdf_high-cdf_low)+d_low;
                }

                beta_n[n]=beta[n];
                Di_n[n]=spreading_function(p, beta[n], wi[n]);

            }

            // p->B133 divided by p->wN must be an interger
            if(p->wN % p->B133 == 0)
                NR=int(p->wN/p->B133);

            // create an angular array of the same length as the frequency array
            count=0;
            for(n = 0; n <NR; ++n)
            {
                for(q = 0; q < p->B133; ++q)
                {
                    beta[count] = beta[q];
                    sinbeta[count] = sin(beta[count]);
                    cosbeta[count] = cos(beta[count]);
                    ++count;
                }
            }

            // randomly re-shuffle the angular array

            if(p->B138>0)
            srand(p->B138_2);

            if(p->B138==0)
            srand((unsigned)time(0));
            
            for(n = 0; n < p->wN; ++n)
            {
                // introduce random index
                int index=rand() % p->wN;
                // keep swaping with the random array to avoid repeating
                double temp=beta[n];
                beta[n] = beta[index];
                beta[index]=temp;
            }

            count=0;
            for(n = 0; n < p->wN; ++n)
            {
                beta[count] = beta[n];
                sinbeta[count] = sin(beta[count]);
                cosbeta[count] = cos(beta[count]);
                ++count;
            }

            count=0;
            for(n = 0; n < p->wN; ++n)
            {
                wi[count] = wi[n];
                dw[count] = dw[n];
                ki[count] = ki[n];
                ++count;
            }

            double* Si_temp;
            p->Darray(Si_temp, p->wN);

            for(n = 0; n < p->wN; ++n)
            Si_temp[n] = Si[n];

            // multiply S with D
            count = 0;
            for(n = 0; n < p->wN; ++n)
            {
                Di[count]=spreading_function(p, beta[count], wi[count]);
                Si[count] = Si_temp[n] * Di[count];
                ++count;
            }

            p->del_Darray(Si_temp, p->wN);
        }
    }
}

double wave_lib_spectrum::spreading_function(lexer* p, double beta, double w)
{
    double D = 0.0;

    // PNJ
    if(p->B130 == 1)
    {
        s_f = p->B134; // CAN BE ANY POSITIVE INTEGER

        if((beta - p->B131) <= PI / 2.0 && (beta - p->B131) >= -PI / 2.0)
            D = (2.0 / PI) * pow(cos(beta - p->B131), s_f);

        if((beta - p->B131) >= PI / 2.0 || (beta - p->B131) <= -PI / 2.0)
            D = 0.0;
    }

    // Mitsuyasu
    if(p->B130 == 2)
    {

        if(p->B134 <= 0.0)  // DON NOT USE THIS OPTION B 135, IT IS WRONG
        {
            if(w <= p->wwp)
                s_f = p->B135 * pow(w / p->wwp, 5.0);

            if(w > p->wwp)
                s_f = p->B135 * pow(w / p->wwp, -2.5);
        }

        if(p->B134 > 0.0)  // USE THIS OPTION B 134, IN THIS CASE, <=85
        {
            s_f = p->B134;
        }

        G_0 = pow(2.0, 2.0 * s_f - 1) / PI * pow(tgamma(s_f + 1.0), 2.0) / tgamma(2.0 * s_f + 1.0);

        D = G_0 * pow(cos((beta - p->B131) / 2.0), 2.0 * s_f);
    }

    return D;
}
