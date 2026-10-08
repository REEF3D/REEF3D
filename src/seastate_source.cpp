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
--------------------------------------------------------------------
Parts of this file are C++ translations of routines of SWAN 41.51:
FAC4WW and SWSNL1 (DIA quadruplets), FAC3WW and SWLTA (LTA triads),
SSURF (Newton linearisation of the Battjes-Janssen dissipation) and
SINTGRL (maximum energy with Battjes-Janssen breaking).

  SWAN (Simulating WAves Nearshore); a third generation wave model
  Copyright (C) 1993-2024  Delft University of Technology

  SWAN is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
--------------------------------------------------------------------*/

#include"seastate_source.h"
#include"seastate_grid.h"
#include"seastate_dispersion.h"
#include<algorithm>
#include<cmath>

namespace
{
    const double pi = 3.14159265358979323846;
    const double tail_p = 4.0;              // spectral tail E ~ sig^-4 (SWAN PWTAIL(1))
    const double rho_ratio = 1.28/1025.0;   // rho_air/rho_water (SWAN PWIND(9))
    const double dia_lambda = 0.25, dia_C = 3.0e7;
    const double dia_cs1 = 5.5, dia_cs2 = 0.833, dia_cs3 = -1.25;
    const double wc_cds = 2.36e-5, wc_stpm = 3.02e-3;
}

seastate_source::seastate_source(const seastate_grid &grid, const seastate_source_param &p) : g(grid), prm(p)
{
    nsig = g.nsig;
    ndir = g.ndir;
    nbin = g.nbin;
    lnr = std::log(g.ratio);

    S.assign(nbin,0.0);
    L.assign(nbin,0.0);
    dNmax.assign(nsig,0.0);
    EL.assign(nsig,0.0);
    TQ.assign(nsig,0.0);
    row.assign(nsig,0.0);

    // powers of the spectral tail (moments, cap, DIA), once
    {
    const double smax = g.sig[nsig-1];
    const double se = smax*std::sqrt(g.ratio);
    pw_smax = std::pow(smax,tail_p);
    pw_se0 = std::pow(se,-tail_p);
    pw_se1 = std::pow(se,1.0-tail_p);
    pw_se2 = std::pow(se,2.0-tail_p);
    pw_fachfr = std::pow(g.ratio,-tail_p);
    }

    // single-frequency grid (surfbeat): no DIA and no triads, so no frequency interpolation
    if(nsig<2)
    return;

    // DIA interpolation (SWAN FAC4WW)
    const double lamm2 = (1.0-dia_lambda)*(1.0-dia_lambda);
    const double lamp2 = (1.0+dia_lambda)*(1.0+dia_lambda);
    const double delth3 = std::acos((lamm2*lamm2 + 4.0 - lamp2*lamp2)/(4.0*lamm2));
    const double delth4 = std::asin(-std::sin(delth3)*lamm2/lamp2);

    dal1 = 1.0/std::pow(1.0+dia_lambda,4.0);
    dal2 = 1.0/std::pow(1.0-dia_lambda,4.0);
    dal3 = 2.0*dal1*dal2;

    const double cidp = std::fabs(delth4/g.dtheta);
    idp  = int(cidp);
    idp1 = idp+1;
    const double widp = cidp - double(idp), widp1 = 1.0-widp;

    const double cidm = std::fabs(delth3/g.dtheta);
    idm  = int(cidm);
    idm1 = idm+1;
    const double widm = cidm - double(idm), widm1 = 1.0-widm;

    const double xis = g.ratio;

    isp  = int(std::log(1.0+dia_lambda)/lnr);
    isp1 = isp+1;
    const double wisp = (1.0+dia_lambda - std::pow(xis,isp))/(std::pow(xis,isp1) - std::pow(xis,isp)), wisp1 = 1.0-wisp;

    ism  = int(std::log(1.0-dia_lambda)/lnr);
    ism1 = ism-1;
    const double wism = (std::pow(xis,ism) - (1.0-dia_lambda))/(std::pow(xis,ism) - std::pow(xis,ism1)), wism1 = 1.0-wism;

    awg[0] = widp *wisp;
    awg[1] = widp1*wisp;
    awg[2] = widp *wisp1;
    awg[3] = widp1*wisp1;
    awg[4] = widm *wism;
    awg[5] = widm1*wism;
    awg[6] = widm *wism1;
    awg[7] = widm1*wism1;

    // E(f,theta) for l = ism1 .. nsig-1-ism1+isp1 (zero below the grid, sig^-4 tail above)
    uoff = -ism1;
    ulen = nsig - 1 - ism1 + isp1 + uoff + 1;
    UE.assign(size_t(ulen)*ndir,0.0);

    // SA1, SA2 for l = -isp1 .. nsig-1-ism1 (zero for l < 0)
    soff = isp1;
    slen = nsig - 1 - ism1 + soff + 1;
    SA1.assign(size_t(slen)*ndir,0.0);
    SA2.assign(size_t(slen)*ndir,0.0);
    DA1C = DA1P = DA1M = DA2C = DA2P = DA2M = SA1;

    af11.assign(nsig-ism1,0.0);
    for(int l=0; l<nsig-ism1; ++l)
    af11[l] = std::pow(g.fmin*std::pow(xis,double(l)),11.0);

    // LTA interpolation at sig/2 (tsm, tsm1) and 2 sig (tsp, tsp1) (SWAN FAC3WW)
    tsm  = int(std::log(0.5)/lnr);
    tsm1 = tsm-1;
    twm  = (std::pow(xis,tsm) - 0.5)/(std::pow(xis,tsm) - std::pow(xis,tsm1));
    twm1 = 1.0-twm;

    tsp  = int(std::log(2.0)/lnr);
    tsp1 = tsp+1;
    twp  = (2.0 - std::pow(xis,tsp))/(std::pow(xis,tsp1) - std::pow(xis,tsp));
    twp1 = 1.0-twp;

    SAL.assign(nsig+tsp1+1,0.0);
}

double seastate_source::ustar_wu(double U10)
{
    const double cd = (U10>7.5) ? (0.8 + 0.065*U10)*1.0e-3 : 1.2875e-3;
    return std::sqrt(cd)*U10;
}

double seastate_source::Qb_bj(double Hrms, double Hm)
{
    const double b = (Hm>0.0 && Hrms>=0.0) ? Hrms/Hm : 0.0;

    if(b<=0.2)
    return 0.0;

    if(b>=1.0)
    return 1.0;

    const double q0 = (b<=0.5) ? 0.0 : (2.0*b-1.0)*(2.0*b-1.0);
    const double b2 = b*b;
    const double z  = std::exp((q0-1.0)/b2);

    return q0 - b2*(q0-z)/(b2-z);
}

void seastate_source::moments(const float *N, double depth, const float *k)
{
    double etot=0.0, actot=0.0, etot1=0.0, edrk=0.0, emax=0.0;

    // row sums: from cap() on the same spectrum, else here
    if(!rows_ready)
    rows(N);
    rows_ready = false;

    for(int l=0; l<nsig; ++l)
    {
    double el = row[l];

    el *= g.sig[l]*g.dtheta;        // E(sig) [m^2 s/rad]

    const double ds = g.dsig[l];
    etot  += el*ds;
    actot += el/g.sig[l]*ds;
    etot1 += el*g.sig[l]*ds;
    edrk  += el/std::sqrt(std::max(double(k[l]),1.0e-12))*ds;

        if(l==nsig-1)
        emax = el;
    }

    // tail E = emax (sig/sigmax)^-4 above the upper edge of the last bin (not for a single-frequency grid)
    if(nsig>1)
    {
    const double smax = g.sig[nsig-1];
    const double se = smax*std::sqrt(g.ratio);
    const double a = emax*pw_smax;
    const double kmax = std::max(double(k[nsig-1]),1.0e-12);

    etot  += a*pw_se1/(tail_p-1.0);
    actot += a*pw_se0/tail_p;
    etot1 += a*pw_se2/(tail_p-2.0);
    edrk  += a*smax/std::sqrt(kmax)*pw_se0/tail_p;
    }

    Etot = etot;
    Hs = 0.0;
    sigm01 = sigm_10 = km_wam = ursell = Qb = 0.0;

    if(etot<=0.0)
    return;

    Hs = 4.0*std::sqrt(etot);
    sigm01 = etot1/etot;
    sigm_10 = etot/actot;
    km_wam = (etot/edrk)*(etot/edrk);

    if(depth>0.0)
    ursell = seastate_gravity*Hs/(2.0*std::sqrt(2.0)*sigm01*sigm01*depth*depth);
}

// sum of n values, eight partial sums (vectorised)
static inline double rowsum(const float *v, int n)
{
    double s[8] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
    int m=0;

    for(; m+8<=n; m+=8)
    for(int k=0; k<8; ++k)
    s[k] += double(v[m+k]);

    for(; m<n; ++m)
    s[0] += double(v[m]);

    return ((s[0]+s[1])+(s[2]+s[3]))+((s[4]+s[5])+(s[6]+s[7]));
}

// sums over the directions per frequency of N (rows), kept for the next moments() of the same spectrum
void seastate_source::rows(const float *N) const
{
    if(g.uniform)
    for(int l=0; l<nsig; ++l)
    row[l] = rowsum(N+g.bin(l,0),ndir);
    else
    for(int l=0; l<nsig; ++l)
    {
    // fine direction sector (A 715): N dth/dtheta
    const float *Nl = N+g.bin(l,0);
    double s = 0.0;
    for(int m=0; m<ndir; ++m)
    s += double(Nl[m])*g.wth[m];
    row[l] = s;
    }
    rows_ready = true;
}

bool seastate_source::cap(float *N, double depth) const
{
    if(!prm.emax || !prm.breaking || prm.breaking_model!=1 || !(depth>0.0))
    return false;

    // total energy with the sig^-4 tail, as moments
    double etot=0.0, elast=0.0;

    if(!rows_ready)
    rows(N);

    for(int l=0; l<nsig; ++l)
    {
    double el = row[l];

    el *= g.sig[l]*g.dtheta;
    etot += el*g.dsig[l];

        if(l==nsig-1)
        elast = el;
    }

    if(nsig>1)
    {
    const double smax = g.sig[nsig-1];
    const double se = smax*std::sqrt(g.ratio);
    etot += elast*pw_smax*pw_se1/(tail_p-1.0);
    }

    const double hm = prm.gamma*depth;
    const double emax = 0.25*hm*hm;

    if(!(etot>emax))
    return false;

    const float f = float(emax/etot);

    for(int b=0; b<nbin; ++b)
    N[b] *= f;

    for(int l=0; l<nsig; ++l)
    row[l] *= double(f);

    return true;
}

void seastate_source::compute(const float *N, double depth, const float *k, const float *cg, double *P, double *D, int ma, int mb, int la, int lb)
{
    moments(N,depth,k);

    wa = ma;
    wb = (mb<0) ? ndir-1 : mb;
    fa = la;
    fb = (lb<0) ? nsig-1 : lb;

    for(int l=fa; l<=fb; ++l)
    for(int m=wa; m<=wb; ++m)
    {
    const int b = g.bin(l,m);
    P[b] = D[b] = S[b] = L[b] = 0.0;
    }

    // DIA writes all bins
    if(prm.dia)
    {
    std::fill(S.begin(),S.end(),0.0);
    std::fill(L.begin(),L.end(),0.0);
    }

    const double grav = seastate_gravity;

    // action density limiter: gamma alpha_PM/(2 sig k^3 c_g)
    if(prm.komen && prm.limiter>0.0)
    for(int l=fa; l<=fb; ++l)
    {
    const double kl = std::max(double(k[l]),1.0e-12);
    dNmax[l] = prm.limiter*0.0081/(2.0*g.sig[l]*kl*kl*kl*std::max(double(cg[l]),1.0e-12));
    }

    // wind input
    ustar = 0.0;

    if(prm.wind && prm.U10>0.0)
    {
    ustar = ustar_wu(prm.U10);

    const double cw = std::cos(prm.wdir), sw = std::sin(prm.wdir);
    const double spm = grav/(28.0*ustar);
    const double ta = prm.Alin/(grav*grav*2.0*pi);

        for(int l=fa; l<=fb; ++l)
        {
        const double sig = g.sig[l];
        const double argu = std::min(2.0,spm/sig);
        const double filter = std::exp(-argu*argu*argu*argu);
        const double freq = sig/(2.0*pi);
        const double reduc = (freq>1.0) ? 1.0/(freq*freq*freq) : 1.0;
        const double cinv = double(k[l])/sig;

            for(int m=wa; m<=wb; ++m)
            {
            const int b = g.bin(l,m);
            const double cosdif = g.costh[m]*cw + g.sinth[m]*sw;

                // linear growth
                if(prm.Alin>0.0 && cosdif>0.0 && sig>=0.7*spm)
                {
                const double t = ustar*cosdif;
                P[b] += reduc*ta/sig*t*t*t*t*filter;
                }

                // exponential growth
                if(prm.komen)
                P[b] += std::max(0.0,0.25*rho_ratio*(28.0*ustar*cinv*cosdif - 1.0))*sig*double(N[b]);
            }
        }
    }

    if(Etot>0.0)
    {
        // whitecapping (Komen)
        if(prm.komen && km_wam>0.0)
        {
        const double stp = km_wam*std::sqrt(Etot)/std::sqrt(wc_stpm);
        const double ck = wc_cds*stp*stp*stp*stp;

            for(int l=fa; l<=fb; ++l)
            {
            const double r = double(k[l])/km_wam;
            const double w = ck*r*sigm_10*r;

            for(int m=wa; m<=wb; ++m)
            D[g.bin(l,m)] += w;
            }
        }

        // depth-induced breaking, Roelvink (1993) for the wave groups of the surfbeat mode (XBeach 'roelvink2'):
        // Qb = 1 - exp(-(H/(gamma h))^n), D/E = 2 alpha f_rep Qb H/h, H = sqrt(8 E) of the instantaneous group
        brk_rate = 0.0;

        if(prm.breaking && prm.breaking_model==2 && depth>0.0)
        {
        const double H = std::sqrt(8.0*Etot);
        const double arg = std::pow(H/(prm.gamma*depth),prm.nroel);

        Qb = std::min(1.0,1.0-std::exp(-std::min(arg,100.0)));
        brk_rate = 2.0*prm.alpha*g.f[0]*Qb*H/depth;

            for(int l=fa; l<=fb; ++l)
            for(int m=wa; m<=wb; ++m)
            D[g.bin(l,m)] += brk_rate;
        }

        // depth-induced breaking (Battjes-Janssen)
        if(prm.breaking && prm.breaking_model==1 && depth>0.0)
        {
        const double Hm = prm.gamma*depth;
        const double bb = 8.0*Etot/(Hm*Hm);

        Qb = Qb_bj(std::sqrt(8.0*Etot),Hm);

        const double ws = (bb<1.0) ? prm.alpha/pi*Qb*sigm01/bb : prm.alpha/pi*sigm01;

        // Newton linearisation of -ws(E) N around the latest spectrum (SWAN SSURF, SbrD):
        // P += sbrd N, D += ws + sbrd with sbrd = ws (1 - Qb)/(bb - Qb) >= 0
        const double sbrd = (bb<1.0 && bb-Qb>1.0e-12) ? ws*(1.0-Qb)/(bb-Qb) : 0.0;

            for(int l=fa; l<=fb; ++l)
            for(int m=wa; m<=wb; ++m)
            {
            const int b = g.bin(l,m);
            D[b] += ws + sbrd;
            P[b] += sbrd*double(N[b]);
            }
        }

        if(prm.dia)
        dia(N,depth);

        if(prm.triads)
        lta(N,depth,k,cg);
    }

    // bottom friction (JONSWAP)
    if(prm.friction && depth>0.0)
    for(int l=fa; l<=fb; ++l)
    {
    const double kd = std::min(30.0,double(k[l])*depth);
    const double s = g.sig[l]/std::sinh(kd);
    const double w = prm.Cb/(grav*grav)*s*s;

        for(int m=wa; m<=wb; ++m)
        D[g.bin(l,m)] += w;
    }

    split(N,P,D);
}

void seastate_source::split(const float *N, double *P, double *D)
{
    for(int l=fa; l<=fb; ++l)
    {
    const int b0 = g.bin(l,wa), nw = wb-wa+1;
    const float *Nl = N+b0;
    const double *Sl = S.data()+b0, *Ll = L.data()+b0;
    double *Pl = P+b0, *Dl = D+b0;

        #pragma GCC ivdep
        for(int n=0; n<nw; ++n)
        {
        const double s = Sl[n], lv = Ll[n], v = double(Nl[n]);
        const double neg = (s<0.0 && v>0.0) ? -s/(v>0.0 ? v : 1.0) : 0.0;
        const double ln = lv<0.0 ? lv : 0.0;

        Pl[n] += (s>0.0 ? s : 0.0) - ln*v;
        Dl[n] += neg - ln;
        }
    }
}

void seastate_source::quadruplets(const float *N, double depth, const float *k, double *Sout, double *dSdN)
{
    moments(N,depth,k);
    wa = 0;
    wb = ndir-1;
    fa = 0;
    fb = nsig-1;
    std::fill(S.begin(),S.end(),0.0);
    std::fill(L.begin(),L.end(),0.0);

    if(Etot>0.0)
    dia(N,depth);

    std::copy(S.begin(),S.end(),Sout);

    if(dSdN!=nullptr)
    std::copy(L.begin(),L.end(),dSdN);
}

void seastate_source::triads(const float *N, double depth, const float *k, const float *cg, double *Sout)
{
    moments(N,depth,k);
    wa = 0;
    wb = ndir-1;
    fa = 0;
    fb = nsig-1;
    std::fill(S.begin(),S.end(),0.0);

    if(Etot>0.0)
    lta(N,depth,k,cg);

    std::copy(S.begin(),S.end(),Sout);
}

// DIA (SWAN SWSNL1): S(f,theta) for E(f,theta) = 2 pi sig N, converted to dN/dt
void seastate_source::dia(const float *N, double depth)
{
    const double grav = seastate_gravity;
    const double fachfr = pw_fachfr;
    const int lhi = nsig - 1 - ism1 + isp1;     // last row of UE
    const int shi = nsig - 1 - ism1;            // last row of SA

    for(int l=ism1; l<=lhi; ++l)
    for(int m=0; m<ndir; ++m)
    {
        if(l<0)
        ue(l,m) = 0.0;
        else if(l<nsig)
        ue(l,m) = 2.0*pi*g.sig[l]*double(N[g.bin(l,m)]);
        else
        ue(l,m) = ue(l-1,m)*fachfr;
    }

    const double x = std::max(0.75*depth*km_wam,0.5);
    const double cons = 1.0/(grav*grav*grav*grav)*(1.0 + dia_cs1/x*(1.0-dia_cs2*x)*std::exp(std::max(-1.0e15,dia_cs3*x)));

    std::fill(SA1.begin(),SA1.end(),0.0);
    std::fill(SA2.begin(),SA2.end(),0.0);
    std::fill(DA1C.begin(),DA1C.end(),0.0);
    std::fill(DA1P.begin(),DA1P.end(),0.0);
    std::fill(DA1M.begin(),DA1M.end(),0.0);
    std::fill(DA2C.begin(),DA2C.end(),0.0);
    std::fill(DA2P.begin(),DA2P.end(),0.0);
    std::fill(DA2M.begin(),DA2M.end(),0.0);

    for(int l=0; l<=shi; ++l)
    for(int m=0; m<ndir; ++m)
    {
    const double e00 = ue(l,m);

        if(e00<=0.0)
        continue;

    const double ep1 = awg[0]*ue(l+isp1,dw(m+idp1)) + awg[1]*ue(l+isp1,dw(m+idp))
                     + awg[2]*ue(l+isp ,dw(m+idp1)) + awg[3]*ue(l+isp ,dw(m+idp));
    const double em1 = awg[4]*ue(l+ism1,dw(m-idm1)) + awg[5]*ue(l+ism1,dw(m-idm))
                     + awg[6]*ue(l+ism ,dw(m-idm1)) + awg[7]*ue(l+ism ,dw(m-idm));
    const double ep2 = awg[0]*ue(l+isp1,dw(m-idp1)) + awg[1]*ue(l+isp1,dw(m-idp))
                     + awg[2]*ue(l+isp ,dw(m-idp1)) + awg[3]*ue(l+isp ,dw(m-idp));
    const double em2 = awg[4]*ue(l+ism1,dw(m+idm1)) + awg[5]*ue(l+ism1,dw(m+idm))
                     + awg[6]*ue(l+ism ,dw(m+idm1)) + awg[7]*ue(l+ism ,dw(m+idm));

    const double factor = cons*af11[l]*e00;
    const double sa1a = e00*(ep1*dal1 + em1*dal2)*dia_C;
    const double sa1b = sa1a - ep1*em1*dal3*dia_C;
    const double sa2a = e00*(ep2*dal1 + em2*dal2)*dia_C;
    const double sa2b = sa2a - ep2*em2*dal3*dia_C;

    sa1(l,m) = factor*sa1b;
    sa2(l,m) = factor*sa2b;

    // derivatives with respect to E00, E+ and E- (SWAN DA1C, DA1P, DA1M, ...)
    const size_t x = sx(l,m);
    DA1C[x] = cons*af11[l]*(sa1a + sa1b);
    DA1P[x] = factor*(dal1*e00 - dal3*em1)*dia_C;
    DA1M[x] = factor*(dal2*e00 - dal3*ep1)*dia_C;
    DA2C[x] = cons*af11[l]*(sa2a + sa2b);
    DA2P[x] = factor*(dal1*e00 - dal3*em2)*dia_C;
    DA2M[x] = factor*(dal2*e00 - dal3*ep2)*dia_C;
    }

    double swg[8];
    for(int n=0; n<8; ++n)
    swg[n] = awg[n]*awg[n];

    // the transfers into the bins of the window of this compute (the quadrant of the sweep)
    for(int l=fa; l<=fb; ++l)
    {
    const double rsigpi = 1.0/(2.0*pi*g.sig[l]);

        for(int m=wa; m<=wb; ++m)
        {
        const double sfnl = -2.0*(sa1(l,m) + sa2(l,m))
            + awg[0]*(sa1(l-isp1,dw(m-idp1)) + sa2(l-isp1,dw(m+idp1)))
            + awg[1]*(sa1(l-isp1,dw(m-idp )) + sa2(l-isp1,dw(m+idp )))
            + awg[2]*(sa1(l-isp ,dw(m-idp1)) + sa2(l-isp ,dw(m+idp1)))
            + awg[3]*(sa1(l-isp ,dw(m-idp )) + sa2(l-isp ,dw(m+idp )))
            + awg[4]*(sa1(l-ism1,dw(m+idm1)) + sa2(l-ism1,dw(m-idm1)))
            + awg[5]*(sa1(l-ism1,dw(m+idm )) + sa2(l-ism1,dw(m-idm )))
            + awg[6]*(sa1(l-ism ,dw(m+idm1)) + sa2(l-ism ,dw(m-idm1)))
            + awg[7]*(sa1(l-ism ,dw(m+idm )) + sa2(l-ism ,dw(m-idm )));

        S[g.bin(l,m)] += sfnl*rsigpi;

        // diagonal derivative dS/dN = d sfnl/dE(f,theta) (SWAN DSNL)
        const double dsnl = -2.0*(DA1C[sx(l,m)] + DA2C[sx(l,m)])
            + swg[0]*(DA1P[sx(l-isp1,dw(m-idp1))] + DA2P[sx(l-isp1,dw(m+idp1))])
            + swg[1]*(DA1P[sx(l-isp1,dw(m-idp ))] + DA2P[sx(l-isp1,dw(m+idp ))])
            + swg[2]*(DA1P[sx(l-isp ,dw(m-idp1))] + DA2P[sx(l-isp ,dw(m+idp1))])
            + swg[3]*(DA1P[sx(l-isp ,dw(m-idp ))] + DA2P[sx(l-isp ,dw(m+idp ))])
            + swg[4]*(DA1M[sx(l-ism1,dw(m+idm1))] + DA2M[sx(l-ism1,dw(m-idm1))])
            + swg[5]*(DA1M[sx(l-ism1,dw(m+idm ))] + DA2M[sx(l-ism1,dw(m-idm ))])
            + swg[6]*(DA1M[sx(l-ism ,dw(m+idm1))] + DA2M[sx(l-ism ,dw(m-idm1))])
            + swg[7]*(DA1M[sx(l-ism ,dw(m+idm ))] + DA2M[sx(l-ism ,dw(m-idm ))]);

        L[g.bin(l,m)] += dsnl;
        }
    }
}

// LTA (original SWAN LTA, SWLTA with ITRIAD 11): self-self sum interactions
void seastate_source::lta(const float *N, double depth, const float *k, const float *cg)
{
    if(!(ursell>=prm.urslim) || depth<=0.0 || sigm01<=0.0)
    return;

    const double grav = seastate_gravity;
    const double biph = 0.5*pi*(std::tanh(prm.urcrit/ursell) - 1.0);
    const double sinbph = std::sin(-biph);

    // last sum frequency below cutfr sig_01
    int ismax = -1;
    for(int l=0; l<nsig; ++l)
    if(g.sig[l]<prm.cutfr*sigm01)
    ismax = l;

    if(ismax<0)
    return;

    // interaction coefficient alpha c cg J^2 (Madsen and Sorensen 1993)
    std::vector<double> &q = TQ;
    std::fill(q.begin(),q.end(),0.0);

    const double d2 = depth*depth, d3 = d2*depth;

    for(int l=0; l<=ismax; ++l)
    if(l+tsm1>=0)
    {
    const double w0  = g.sig[l];
    const double wn0 = double(k[l]);
    const double c0  = w0/wn0;
    const double wm  = twm*g.sig[l+tsm1] + twm1*g.sig[l+tsm];
    const double wnm = twm*double(k[l+tsm1]) + twm1*double(k[l+tsm]);

        if(wnm<=0.0 || wn0<=0.0)
        continue;

    const double R  = (0.5 + wm*wm/(grav*depth*wnm*wnm))*(2.0*wnm)*(2.0*wnm);
    const double Sd = (-2.0/grav)*(grav*depth*wn0 + 2.0/15.0*grav*d3*wn0*wn0*wn0 - 2.0/5.0*d2*wn0*w0*w0);

        if(Sd==0.0)
        continue;

    const double J = R/Sd;
    q[l] = prm.alphaEB*double(cg[l])*c0*J*J*sinbph;
    }

    for(int m=wa; m<=wb; ++m)
    {
        for(int l=0; l<nsig; ++l)
        EL[l] = 2.0*pi*g.sig[l]*double(N[g.bin(l,m)]);

        std::fill(SAL.begin(),SAL.end(),0.0);

        for(int l=0; l<=ismax; ++l)
        if(l+tsm1>=0)
        {
        const double e0 = EL[l];
        const double em = twm*EL[l+tsm1] + twm1*EL[l+tsm];

        SAL[l] = std::max(0.0,q[l]*(em*(em-e0) - e0*em));
        }

        for(int l=fa; l<=fb; ++l)
        {
        const double stri = SAL[l] - 2.0*(twp*SAL[l+tsp1] + twp1*SAL[l+tsp]);

        S[g.bin(l,m)] += stri/(2.0*pi*g.sig[l]);
        }
    }
}
