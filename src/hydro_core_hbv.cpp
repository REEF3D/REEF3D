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

#include"hydro_core.h"
#include"lexer.h"
#include<algorithm>
#include<cmath>

/*
HBV per owned cell for one step dt (depths in mm over the area they refer to, dtd = dt/1 day):

band b (elevation zb, area fraction bf, glacier fraction g), station forcing corrected to zb:
    T = sum_s w_s (T_s - TCALT/100 (zb - z_s)),  Pr = sum_s w_s P_s max(0, 1 + PCALT/10000 (zb - z_s)) dt
    PET from the columns (weights as P without the elevation correction) or Oudin with T
snow:   T <  TT: snowfall SF = Pr PCORR SFCF, rain 0;  T >= TT: rain = Pr PCORR
        melt = min(CFMAX dtd (T-TT), SP) for T > TT,  refreezing = min(CFR CFMAX dtd (TT-T), WC) for T < TT
        SP += SF + refreezing - melt;  WC += rain + melt - refreezing;  out = max(WC - CWH SP, 0);  WC -= out
glacier (fraction g of the band, snow free after the snow routine):
        icemelt = CFGLAC CFMAX dtd max(T-TT, 0), out and icemelt go to the upper zone directly
soil (fraction 1-g):
        r = out (SM/FC)^BETA;  SM += out - r;  SM > FC: r += SM - FC, SM = FC
        Ea = 0 under snow (SP > 0), otherwise min(PET min(SM/(LP FC), 1), SM);  SM -= Ea
response per cell, recharge R = sum_b bf [(1-g) r + g (out + icemelt)]:
        SUZ += R;  perc = min(PERC dtd, SUZ);  SUZ -= perc;  SLZ += perc
        Q0 = k0 max(SUZ - UZL, 0), SUZ -= Q0;  Q1 = k1 SUZ, SUZ -= Q1;  Q2 = k2 SLZ, SLZ -= Q2
        with k = 1 - (1-K)^dtd
MAXBAS: Q0+Q1+Q2 is spread over the following steps with the triangular weights
*/
void hydro_core::vertical(lexer *p)
{
    const int ns=int(sx.size());
    const double dtd=dt/86400.0;
    const double dth=dt/3600.0;
    const double tday=(t+0.5*dt)/86400.0;

    double vp=0.0, vea=0.0, vgm=0.0;

    for(int n=n0; n<n1; ++n)
    {
        const int o=n-n0;
        const double *w=&wgt[size_t(o)*ns];

        const double tt=P(n,TT), cfmax=P(n,CFMAX)*dtd, fc=P(n,FC), lp=P(n,LP), beta=P(n,BETA);
        double recharge=0.0, sP=0.0, sT=0.0, sE=0.0, sEa=0.0, sgm=0.0;

        for(int b=0; b<nb; ++b)
        {
            const int q=o*nb+b;
            const double zb=bz[n*nb+b], fb=bf[n*nb+b], g=bg[n*nb+b];
            if(fb<=0.0)
            continue;

            double T=0.0, Pr=0.0, E=0.0;
            for(int k=0; k<ns; ++k)
            {
                T  += w[k]*(mT[k] - P(n,TCALT)/100.0*(zb-sz[k]));
                Pr += w[k]*mP[k]*max(0.0,1.0+P(n,PCALT)/10000.0*(zb-sz[k]));
                E  += w[k]*mE[k];
            }
            Pr*=dth;
            E = fhasE ? E*dth : pet_oudin(p,T,tday)*dtd;

            // snow
            double sf=0.0, rain=0.0;
            if(T<tt)
            sf=Pr*P(n,PCORR)*P(n,SFCF);
            else
            rain=Pr*P(n,PCORR);

            const double melt = (T>tt) ? min(cfmax*(T-tt),SP[q]) : 0.0;
            const double refr = (T<tt) ? min(P(n,CFR)*cfmax*(tt-T),WC[q]) : 0.0;
            SP[q] += sf + refr - melt;
            WC[q] += rain + melt - refr;
            const double out = max(WC[q]-P(n,CWH)*SP[q],0.0);
            WC[q] -= out;

            // glacier
            const double icemelt = (g>0.0 && SP[q]<=0.0) ? P(n,CFGLAC)*cfmax*max(T-tt,0.0) : 0.0;

            // soil
            double r=0.0, ea=0.0;
            if(g<1.0)
            {
                r = out*pow(min(SM[q]/fc,1.0),beta);
                SM[q] += out - r;
                if(SM[q]>fc)
                {
                    r += SM[q]-fc;
                    SM[q] = fc;
                }

                if(SP[q]<=0.0)
                {
                    ea = min(E*min(SM[q]/(lp*fc),1.0),SM[q]);
                    SM[q] -= ea;
                }
            }

            recharge += fb*((1.0-g)*r + g*(out+icemelt));
            sP  += fb*(sf+rain);
            sT  += fb*T;
            sE  += fb*E;
            sEa += fb*(1.0-g)*ea;
            sgm += fb*g*icemelt;
        }

        // response
        double &suz=SUZ[o], &slz=SLZ[o];
        suz += recharge;
        const double perc=min(P(n,PERC)*dtd,suz);
        suz -= perc;
        slz += perc;

        const double k0=1.0-pow(1.0-P(n,K0),dtd), k1=1.0-pow(1.0-P(n,K1),dtd), k2=1.0-pow(1.0-P(n,K2),dtd);
        const double q0=k0*max(suz-P(n,UZL),0.0);
        suz -= q0;
        const double q1=k1*suz;
        suz -= q1;
        const double q2=k2*slz;
        slz -= q2;

        // MAXBAS
        double *bb=&buf[size_t(o)*nmb];
        const double *wm=&mbw[size_t(n)*nmb];
        const double qs=q0+q1+q2;
        for(int k=0; k<nmb; ++k)
        bb[k]+=wm[k]*qs;
        const double qout=bb[0];
        for(int k=0; k<nmb-1; ++k)
        bb[k]=bb[k+1];
        bb[nmb-1]=0.0;

        qgen[n] = qout*1.0e-3*area[n]/dt;

        cP[o]=sP/dth;
        cT[o]=sT;
        cE[o]=sE/dth;
        cEa[o]=sEa/dth;
        cRun[o]=qout/dth;

        vp  += sP*1.0e-3*area[n];
        vea += sEa*1.0e-3*area[n];
        vgm += sgm*1.0e-3*area[n];
    }

    Vp+=vp;
    Vea+=vea;
    Vgm+=vgm;
}

// water stored in the owned cells [m3]: snow, soil, upper and lower zone, MAXBAS buffer
double hydro_core::owned_storage() const
{
    double s=0.0;
    for(int n=n0; n<n1; ++n)
    {
        const int o=n-n0;
        double d=SUZ[o]+SLZ[o];
        for(int b=0; b<nb; ++b)
        {
            const int q=o*nb+b;
            d += bf[n*nb+b]*(SP[q]+WC[q]+(1.0-bg[n*nb+b])*SM[q]);
        }
        for(int k=0; k<nmb; ++k)
        d += buf[size_t(o)*nmb+k];

        s += d*1.0e-3*area[n];
    }
    return s;
}
