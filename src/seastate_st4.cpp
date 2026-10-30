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
This file contains a C++ translation of routines of the ECMWF wave
model ecWAM 1.5.14 (https://github.com/ecmwf-ifs/ecwam):

  (C) Copyright 1989- ECMWF.

  This software is licensed under the terms of the Apache Licence
  Version 2.0 which can be obtained at
  http://www.apache.org/licenses/LICENSE-2.0.
  In applying this licence, ECMWF does not waive the privileges and
  immunities granted to it by virtue of its status as an
  intergovernmental organisation nor does it submit to any
  jurisdiction.

Translated: SINPUT_ARD (NGST 1, no LLNORMAGAM), SDISSIP_ARD (SSDSC3 0,
SSDSC5 0), TAUT_Z0 (without LLGCBZ0), STRESSO and TAU_PHI_HF (without
LLGCBZ0 and LLNORMAGAM), FRCUTINDEX (no ice), FEMEANWS, FKMEAN (FM1),
CHNKMIN, CDM, INIT_X0TAUHF, INIT_SDISS_ARDH, TABU_SWELLFT, KERKEI,
KZEONE (series branch, |z| < 4.3), IMPHFTAIL, parameters of SETWAVPHYS
(IPHYS 1, LLGCBZ0 and LLNORMAGAM false).
Changes: one grid point; the SEASTATE spectral grid (its bin widths
DFIM = dsig dtheta / 2 pi, directions from 0, the compact direction
sector without wrapping); action density N(sig,theta) with F(f,theta)
= 2 pi sig N; Cartesian directions (x = cos, y = sin); no gustiness
(WSTAR 0); the lag of the wave stress by a time step replaced by the
fixed point of u* and tau_w for the spectrum of the cell; the source
terms returned as P = gamma N and D >= 0 for the implicit solver; the
diagnostic tail imposed after each solve of a cell; gravity of
SEASTATE (9.81).
--------------------------------------------------------------------*/

#include"seastate_st4.h"
#include"seastate_grid.h"
#include<algorithm>
#include<cmath>

namespace
{
    const double pi = 3.14159265358979323846;
    const double zpi = 2.0*pi;
    const double xkappa = 0.40, xnlev = 10.0;
    const double epsmin = 0.1e-32, epsus = 1.0e-6, eps1 = 0.00001;
    const double acd = 8.0e-4, bcd = 8.0e-5, cdmax = 0.0025;
    const double swellf = 0.66, swellf2 = -0.018, swellf3 = 0.022, swellf5 = 1.2, swellf6 = 1.0;
    const double abmin = 0.3, abmax = 8.0;
    const int iab = 200;
    const int isdsdth = 80;
    const double fric = 28.0, wetail = 0.25, frtail = 0.2;
    const int jtot = 19;

    // Fortran NINT
    inline int nint(double x) {return x>=0.0 ? int(std::floor(x+0.5)) : -int(std::floor(-x+0.5));}
}

// exp(z) K0(z), exp(z) K1(z), z = x + i y (KZEONE, series for |z|^2 < 18.49)
void seastate_st4::kzeone(double x, double y, double &re0, double &im0, double &re1, double &im1)
{
    double r2 = x*x + y*y;
    double x2 = x/2.0, y2 = y/2.0;
    double p1 = x2*x2, p2 = y2*y2;
    double t1 = -(std::log(p1+p2)/2.0 + 0.5772156649015329);
    const double t2 = -std::atan2(y,x);
    x2 = p1 - p2;
    y2 = x*y2;
    double rterm = 1.0, iterm = 0.0;
    re0 = t1;
    im0 = t2;
    t1 = t1 + 0.5;
    re1 = t1;
    im1 = t2;
    p2 = std::sqrt(r2);
    double L = 2.106*p2 + 4.4;
    if(p2<8.0e-1)
    L = 2.129*p2 + 4.0;
    const int ll = nint(L);

    for(int n=1; n<=ll; ++n)
    {
    p1 = double(n);
    p2 = double(n*n);
    const double r1 = rterm;
    rterm = (r1*x2 - iterm*y2)/p2;
    iterm = (r1*y2 + iterm*x2)/p2;
    t1 = t1 + 0.5/p1;
    re0 = re0 + t1*rterm - t2*iterm;
    im0 = im0 + t1*iterm + t2*rterm;
    p1 = p1 + 1.0;
    t1 = t1 + 0.5/p1;
    re1 = re1 + (t1*rterm - t2*iterm)/p1;
    im1 = im1 + (t1*iterm + t2*rterm)/p1;
    }

    const double r1 = x/r2 - 0.5*(x*re1 - y*im1);
    r2 = -y/r2 - 0.5*(x*im1 + y*re1);
    p1 = std::exp(x);
    re0 = p1*re0;
    im0 = p1*im0;
    re1 = p1*r1;
    im1 = p1*r2;
}

void seastate_st4::kerkei(double x, double &ker, double &kei)
{
    const double zr = x*0.5*std::sqrt(2.0);
    double cyr, cyi, c1, c2;
    kzeone(zr,zr,cyr,cyi,c1,c2);
    ker = cyr/std::exp(zr);
    kei = cyi/std::exp(zr);
}

double seastate_st4::chnkmin(double u10) const
{
    return prm.alphamin + (prm.alpha-prm.alphamin)*0.5*(1.0 - std::tanh(u10-prm.chnkmin_u));
}

void seastate_st4::init(const seastate_grid &g, const seastate_st4_param &p, double gravity)
{
    pg = &g;
    prm = p;
    grav = gravity;
    nsig = g.nsig;
    ndir = g.ndir;
    nbin = g.nbin;
    rnum = 0.11*prm.nu;
    raorw = std::max(prm.rhoa,1.0)/prm.rhow;

    gam.assign(nbin,0.0);
    dst.assign(nbin,0.0);
    ustl.assign(nsig,0.0);
    udl.assign(nsig,0.0);
    sxl.assign(nsig,0.0);
    syl.assign(nsig,0.0);
    wsl.assign(nsig,0.0);
    bth.assign(ndir,0.0);
    wcl.assign(ndir,0.0);

    // DFIM and RHOWG_DFIM on the SEASTATE grid
    const double dth = g.dtheta;
    const double xlr = std::log(g.ratio);
    dfim.assign(nsig,0.0);
    rhowg.assign(nsig,0.0);
    for(int l=0; l<nsig; ++l)
    {
    dfim[l] = g.dsig[l]*dth/zpi;
    rhowg[l] = prm.rhow*grav*dth*xlr*g.f[l];
    }
    rhowg[0] *= 0.5;
    rhowg[nsig-1] *= 0.5;

    // FLMAX = (ALPHAPMAX/pi)/((2 pi)^4/g^2 f^5) of F(f,theta), as N = F/(2 pi sig)
    nmax.assign(nsig,0.0);
    for(int l=0; l<nsig; ++l)
    nmax[l] = (prm.alphapmax/pi)/(std::pow(zpi,4.0)/(grav*grav)*std::pow(g.f[l],5.0))/(zpi*g.sig[l]);

    // INIT_X0TAUHF (LLCAPCHNK: ALPHAMIN)
    {
    const double bk2 = prm.betamax/(xkappa*xkappa);
    const double alph = prm.capchnk ? prm.alphamin : prm.alpha;
    double x0 = 0.005;
        for(int j=0; j<30; ++j)
        {
        const double ff = std::exp(xkappa/(x0+prm.zalp));
        const double f = alph*x0*x0*ff - 1.0;
        if(f==0.0)
        break;
        const double df = alph*ff*(2.0*x0 - xkappa*(x0/(x0+prm.zalp))*(x0/(x0+prm.zalp)));
        x0 = x0 - f/df;
        }
    x0tauhf = x0;
    const double c1 = bk2/3.0;
    wtauhf[0] = c1;
    for(int j=2; j<=jtot-1; j+=2)
    {
    wtauhf[j-1] = 4.0*c1;
    wtauhf[j] = 2.0*c1;
    }
    wtauhf[jtot-1] = c1;
    }

    // INIT_SDISS_ARDH: +-80 deg, cos^2, the end bins truncated
    {
    const double rad = pi/180.0;
    const int nangd = ndir/2;
    nsat = std::min(nint(isdsdth*rad/dth),nangd-1);
    double dtr = isdsdth*rad - (double(nsat)*dth - 0.5*dth);
    dtr = std::max(0.0,std::min(dtr,dth));
    const int nw = 2*nsat+1;
    satidx.assign(size_t(ndir)*nw,-1);
    satw.assign(size_t(ndir)*nw,0.0);
        for(int k=0; k<ndir; ++k)
        for(int i=k-nsat; i<=k+nsat; ++i)
        {
        int j = i;
            if(g.periodic)
            {
            if(i<0) j = i+ndir;
            if(i>=ndir) j = i-ndir;
            }
            else if(i<0 || i>=ndir)
            continue;
        const double dl = (i==k-nsat || i==k+nsat) ? dtr : dth;
        const double c = std::cos(g.theta[k]-g.theta[j]);
        satidx[size_t(k)*nw + (i-(k-nsat))] = j;
        satw[size_t(k)*nw + (i-(k-nsat))] = dl*c*c;
        }
    }

    // TABU_SWELLFT
    swellft.assign(iab+1,0.0);
    {
    double dzeta0 = 0.0;
    const double delab = (abmax-abmin)/double(iab);
    const double l10 = std::log(10.0);
        for(int i=1; i<=iab; ++i)
        {
        const double abrlog = abmin + double(i)*delab;
        const double abr = std::exp(abrlog*l10);
        const double fact = 1.0/abr/(21.2*xkappa);
        double fsubw = 0.05;
            for(int it=0; it<100; ++it)
            {
            const double fsubwmemo = fsubw, dzeta0memo = dzeta0;
            dzeta0 = fact*std::pow(fsubw,-0.5);
            double ker, kei;
            kerkei(2.0*std::sqrt(dzeta0),ker,kei);
            fsubw = 0.08/(ker*ker + kei*kei);
            fsubw = 0.5*(fsubwmemo + fsubw);
            dzeta0 = 0.5*(dzeta0memo + dzeta0);
            }
        swellft[i] = fsubw;
        }
    }
}

// TAUT_Z0 without the gravity-capillary model
void seastate_st4::taut_z0(int iusfg, double utop, double udir, double tw, double twdir, double &ustar, double &zz0, double &zz0b, double &chr) const
{
    const int niter = 18;
    const double twoxmp1 = 3.0;
    const double xlogxl = std::log(xnlev);
    const double us2totauw = 1.0 + eps1;
    const double gm1 = 1.0/grav;

    const double cosdiff = std::cos(udir-twdir);
    const double tauwact = std::max(tw*cosdiff,epsmin);
    const double tauweff = tauwact*us2totauw;

    double xmin, alphaog;
    if(prm.capchnk)
    {
    const double cm = chnkmin(utop);
    xmin = 0.15*(prm.alpha-cm);
    alphaog = cm*gm1;
    }
    else
    {
    xmin = 0.0;
    alphaog = prm.alpha*gm1;
    }

    const double xkutop = xkappa*utop;
    const double ustold = (1-iusfg)*utop*std::sqrt(std::min(acd+bcd*utop,cdmax)) + iusfg*ustar;
    double tauold = std::max(ustold*ustold,tauweff);
    ustar = std::sqrt(tauold);
    double ustm1 = 1.0/std::max(ustar,epsus);
    double z0ch = 0.0;

    for(int it=0; it<niter; ++it)
    {
    const double x = std::max(tauwact/tauold,xmin);
    z0ch = alphaog*tauold/std::sqrt(1.0-x);
    const double z0vis = rnum*ustm1;
    const double z0tot = z0ch + z0vis;
    const double xologz0 = 1.0/(xlogxl - std::log(z0tot));
    const double f = ustar - xkutop*xologz0;
    const double zz = ustm1*(z0ch*(2.0-twoxmp1*x)/(1.0-x) - z0vis)/z0tot;
    const double delf = 1.0 - xkutop*xologz0*xologz0*zz;
    if(delf!=0.0)
    ustar = ustar - f/delf;
    const double taunew = std::max(ustar*ustar,tauweff);
    ustar = std::sqrt(taunew);
    if(taunew==tauold)
    break;
    ustm1 = 1.0/std::max(ustar,epsus);
    tauold = taunew;
    }

    zz0 = z0ch;
    zz0b = alphaog*tauold;
    chr = std::max(grav*zz0*ustm1*ustm1,prm.alphamin);
}

// SINPUT_ARD for all bins: growth rates gam, swell damping dst (sneg), the stress of the input per frequency
// (sxl, syl: sum of gamma F cos, sin) and the wind-sea sums (wsl)
void seastate_st4::input(const float *N, const float *k, const float *cg, double us, double zz0, double wdir, bool sneg)
{
    const seastate_grid &g = *pg;
    const double cw = std::cos(wdir), sw = std::sin(wdir);
    const double shelter = std::fabs(prm.tauwshelter);
    const bool lshelter = (shelter!=0.0);
    const double const1 = prm.betamax/(xkappa*xkappa);
    const double rogoroair = grav/raorw;
    const double fu = std::fabs(swellf3), fud = swellf2;
    const double nu_air = prm.nu;

    double taux = us*us*cw, tauy = us*us*sw, xs = 0.0, ys = 0.0;
    double ustp = us, udir = wdir;
    const double aird_pvisc = sw_pvisc*raorw;

    for(int l=0; l<nsig; ++l)
    {
    const double sig = g.sig[l], kl = double(k[l]);
    const double cinv = kl/sig;
    double constf = 0.0;

        if(lshelter)
        {
        const double taupx = taux - shelter*xs, taupy = tauy - shelter*ys;
        udir = std::atan2(taupy,taupx);
        ustp = std::pow(taupx*taupx + taupy*taupy,0.25);
        constf = rogoroair*cinv*dfim[l];
        }

    ustl[l] = ustp;
    udl[l] = udir;

    const double ucn = ustp*cinv;
    const double ucnzalpd = xkappa/(ucn + prm.zalp);
    const double zcn = std::log(kl*zz0);
    const double cnsn = sig*const1*raorw;
    const double cl = std::cos(udir), sl = std::sin(udir);
    const double fs = zpi*sig;

    double coef = 0.0, dstab1 = 0.0, temp1 = 0.0;
        if(sneg)
        {
        coef = -swellf*16.0*sig*sig/grav;
        const double coef5 = -swellf5*2.0*std::sqrt(2.0*nu_air*sig);
        dstab1 = coef5*aird_pvisc*kl;
        temp1 = coef*raorw;
        }

    double sx = 0.0, sy = 0.0, ws = 0.0, xa = 0.0, ya = 0.0;
    const float *Nl = N + g.bin(l,0);
    double *gl = gam.data() + g.bin(l,0), *dl = dst.data() + g.bin(l,0);

        for(int m=0; m<ndir; ++m)
        {
        const double coslp = g.costh[m]*cl + g.sinth[m]*sl;
        const double F = fs*double(Nl[m]);
        double gm = 0.0;

            if(coslp>0.01)
            {
            const double x = coslp*ucn;
            const double zlog = zcn + ucnzalpd/coslp;
                if(zlog<0.0)
                {
                const double zlog2x = zlog*zlog*x;
                gm = std::exp(zlog)*zlog2x*zlog2x*cnsn;
                ws += F;
                }
            }

        gl[m] = gm;
        const double slp = gm*F;
        sx += slp*g.costh[m];
        sy += slp*g.sinth[m];

            if(sneg)
            {
            const double dstab2 = temp1*(sw_temp2 + (fu + fud*coslp)*ustp);
            dl[m] = -(dstab1 + sw_pturb*dstab2);
            }
            else
            dl[m] = 0.0;
        }

    sxl[l] = sx;
    syl[l] = sy;
    wsl[l] = ws;

        if(lshelter)
        {
        xa = sx*constf;
        ya = sy*constf;
        xs += xa;
        ys += ya;
        }
    }
}

// FKMEAN (FM1) and FEMEANWS with the tails of ecWAM, FRCUTINDEX
void seastate_st4::cutoff(const float *N, double us)
{
    const seastate_grid &g = *pg;
    double em = epsmin, fmv = epsmin, ew = epsmin, fw = epsmin, temp = 0.0, tw = 0.0;

    for(int l=0; l<nsig; ++l)
    {
    // sum of F over the directions
    double s = 0.0;
    const float *Nl = N + g.bin(l,0);
    for(int m=0; m<ndir; ++m)
    s += double(Nl[m]);
    temp = zpi*g.sig[l]*s;
    tw = wsl[l];
    em += dfim[l]*temp;
    fmv += dfim[l]/g.f[l]*temp;
    ew += dfim[l]*tw;
    fw += dfim[l]/g.f[l]*tw;
    }

    const double delt25 = wetail*g.f[nsig-1]*g.dtheta, delt2 = frtail*g.dtheta;
    em += delt25*temp;
    fmv += delt2*temp;
    ew += delt25*tw;
    fw += delt2*tw;
    fm = em/fmv;
    fmws = ew/fw;

    const double fpmh = prm.tailfactor/g.f[0];
    const double fppm = prm.tailfactor_pm*grav/(fric*g.sig[0]);
    const double fm2 = std::max(fmws,fm)*fpmh;
    const double fpm = fppm/std::max(us,epsmin);
    const double fpm4 = std::max(fm2,fpm);
    int mij = nint(std::log10(fpm4)/std::log10(g.ratio)) + 1;
    mij = std::min(std::max(1,mij),nsig);
    lc = mij-1;
}

// TAU_PHI_HF: stress of the f^-5 tail above the cut-off m (no gravity-capillary model, no renormalisation)
double seastate_st4::tau_hf(int m, bool shelter, double us, double zz0, double wdir, const float *N, double ust) const
{
    const seastate_grid &g = *pg;
    const double gm1 = 1.0/grav;
    const double zsupmax = 0.0;
    const double x0g = x0tauhf*grav;
    const double xloggz0 = std::log(grav*zz0);
    const double omegacc = std::max(g.sig[m],x0g/ust);
    const double sqrtz0og = std::sqrt(zz0*gm1);
    const double sqrtgz0 = 1.0/sqrtz0og;
    const double yc = omegacc*sqrtz0og;
    const double zinf = std::log(yc);
    const double fr5 = std::pow(g.f[m],5.0);
    const double consttau = std::pow(zpi,4.0)*gm1*gm1*fr5;
    const double cw = std::cos(wdir), sw = std::sin(wdir);

    double f1dcos3 = 0.0;
    const float *Nm = N + g.bin(m,0);
    for(int k=0; k<ndir; ++k)
    {
    const double cosw = std::max(g.costh[k]*cw + g.sinth[k]*sw,0.0);
    f1dcos3 += zpi*g.sig[m]*double(Nm[k])*cosw*cosw*cosw;
    }
    f1dcos3 *= g.dtheta;

    double taul = ust*ust;
    const double delz = std::max((zsupmax-zinf)/double(jtot-1),0.0);
    double tauhf = 0.0;

    if(shelter)
    {
        for(int j=0; j<jtot; ++j)
        {
        const double y = std::exp(zinf + double(j)*delz);
        const double omega = y*sqrtgz0;
        const double cm1 = omega*gm1;
        const double zx = ust*cm1 + prm.zalp;
        const double zarg = xkappa/zx;
        double zlog = xloggz0 + 2.0*std::log(cm1) + zarg;
        zlog = std::min(zlog,0.0);
        const double zbeta = zlog*zlog*zlog*zlog*std::exp(zlog);
        const double fnc2 = f1dcos3*consttau*zbeta*taul*wtauhf[j]*delz;
        taul = std::max(taul - prm.tauwshelter*fnc2,0.0);
        ust = std::sqrt(taul);
        tauhf += fnc2;
        }
    }
    else
    {
        for(int j=0; j<jtot; ++j)
        {
        const double y = std::exp(zinf + double(j)*delz);
        const double omega = y*sqrtgz0;
        const double cm1 = omega*gm1;
        const double zx = ust*cm1 + prm.zalp;
        const double zarg = xkappa/zx;
        double zlog = xloggz0 + 2.0*std::log(cm1) + zarg;
        zlog = std::min(zlog,0.0);
        const double zbeta = zlog*zlog*zlog*zlog*std::exp(zlog);
        tauhf += zbeta*wtauhf[j];
        }
    tauhf = f1dcos3*consttau*taul*tauhf*delz;
    }

    (void)us;
    return tauhf;
}

// STRESSO: the wave stress of the prognostic part (up to the cut-off lc, half there) and of the tail
void seastate_st4::stress(const float *N, const float *k, double us, double zz0, double wdir, double &tw, double &twdir) const
{
    const seastate_grid &g = *pg;
    double xs = 0.0, ys = 0.0;

    for(int l=0; l<=lc; ++l)
    {
    double w = rhowg[l];
    if(l==lc && lc!=nsig-1)
    w *= 0.5;
    const double c = w*double(k[l])/g.sig[l];
    xs += c*sxl[l];
    ys += c*syl[l];
    }

    xs /= std::max(prm.rhoa,1.0);
    ys /= std::max(prm.rhoa,1.0);

    const bool shelter = (prm.tauwshelter!=0.0);
    double usdirp = wdir, ust = us;

    if(shelter)
    {
    const double taupx = us*us*std::cos(wdir) - prm.tauwshelter*xs;
    const double taupy = us*us*std::sin(wdir) - prm.tauwshelter*ys;
    usdirp = std::atan2(taupy,taupx);
    ust = std::pow(taupx*taupx + taupy*taupy,0.25);
    }

    const double tauhf = tau_hf(lc,shelter,us,zz0,wdir,N,ust);
    xs += tauhf*std::cos(usdirp);
    ys += tauhf*std::sin(usdirp);
    tw = std::max(std::sqrt(xs*xs + ys*ys),0.0);
    twdir = std::atan2(ys,xs);
    tw = std::min(tw,us*us/(1.0+eps1));
}

void seastate_st4::wind(const float *N, const double *row, const float *k, const float *cg, double U10, double wdir, bool freeze)
{
    if(freeze)
    return;

    calmed = false;
    const seastate_grid &g = *pg;

    // the fixed point of u* and tau_w (ecWAM: tau_w of the last time step, two calls of SINFLX): tau_w = H(tau_w),
    // H the wave stress of the input with u*, z0 of TAUT_Z0(tau_w); secant steps on H(tau_w) - tau_w (the plain
    // iteration where the secant leaves [0, u*^2)), from tau_w = 0, until u* changes by less than tol
    double tw = 0.0, twd = wdir, us = 0.0, prev = 0.0;
    double zz0 = 0.0, zz0b = 0.0, chr = 0.0;
    double t0 = 0.0, h0 = 0.0;
    int it = 0;

    for(it=0; it<prm.maxit; ++it)
    {
    taut_z0(it>0 ? 1 : 0,U10,wdir,tw,twd,us,zz0,zz0b,chr);
    if(it>0 && std::fabs(us-prev)<=prm.tol*us)
    break;
    prev = us;

    double h = 0.0, hd = wdir;
    input(N,k,cg,us,zz0,wdir,false);
    cutoff(N,us);
    stress(N,k,us,zz0,wdir,h,hd);

    double tn = h;
        if(it>0)
        {
        const double f1 = h-tw, f0 = h0-t0;
            if(std::fabs(f1-f0)>1.0e-14*std::max(us*us,1.0e-12))
            {
            const double ts = tw - f1*(tw-t0)/(f1-f0);
            if(ts>0.0 && ts<us*us)
            tn = ts;
            }
        }
    t0 = tw;
    h0 = h;
    tw = tn;
    twd = hd;
    }
    iters = it;
    if(it==prm.maxit)
    taut_z0(1,U10,wdir,tw,twd,us,zz0,zz0b,chr);

    // swell damping of the cell (SINPUT_ARD with LLSNEG): orbital velocity, amplitude, Reynolds number
    if(prm.swell)
    {
    double uorbt = epsmin, aorb = epsmin;
        for(int l=0; l<nsig; ++l)
        {
        const double sig = g.sig[l];
        const double temp = zpi*sig*row[l];      // sum of F over the directions
        uorbt += dfim[l]*sig*sig*temp;
        aorb += dfim[l]*temp;
        }
    uorbt = 2.0*std::sqrt(uorbt);
    aorb = 2.0*std::sqrt(aorb);
    const double re = 4.0/prm.nu*uorbt*aorb;
    const double z0vis = rnum/std::max(us,0.0001);
    const double z0tub = prm.z0rat*std::min(prm.z0tubmax,zz0);
    const double z0noz = std::max(z0vis,z0tub);
    const double zorb = aorb/z0noz;
    const double delabm1 = double(iab)/(abmax-abmin);
    const double xi = (std::log10(std::max(zorb,3.0)) - abmin)*delabm1;
    const int ind = std::min(iab-1,int(xi));
    const double deli1 = std::min(1.0,xi-double(ind)), deli2 = 1.0-deli1;
    const double fww = swellft[ind]*deli2 + swellft[ind+1]*deli1;
    sw_temp2 = fww*uorbt;
    const double re_c = (swellf6==1.0) ? prm.swellf4 : prm.swellf4*std::pow(2.0/aorb,1.0-swellf6);
        if(prm.swellf7>0.0)
        {
        const double smooth = 0.5*std::tanh((re-re_c)/prm.swellf7);
        sw_pturb = 0.5 + smooth;
        sw_pvisc = 0.5 - smooth;
        }
        else
        {
        sw_pturb = (re<=re_c) ? 0.0 : 0.5;
        sw_pvisc = (re<=re_c) ? 0.5 : 0.0;
        }
    }

    input(N,k,cg,us,zz0,wdir,prm.swell);
    cutoff(N,us);
    stress(N,k,us,zz0,wdir,tw,twd);

    ust = us;
    z0 = zz0;
    z0b = zz0b;
    chrnck = chr;
    tauw = tw;
    tauwdir = twd;
}

// SDISSIP_ARD of frequency l: whitecapping rate (>= 0) per direction
void seastate_st4::dissipation(const float *N, const float *k, const float *cg, int l, double *Dl)
{
    const seastate_grid &g = *pg;
    const double kl = double(k[l]);
    const double sig = g.sig[l];
    const double facsat = kl*kl*kl*double(cg[l])*sig;      // k^3 cg E(sig) = k/(2 pi) k^2 cg F(f)
    const double tmp03 = 1.0/prm.sdsbr;
    const int nw = 2*nsat+1;
    const float *Nl = N + g.bin(l,0);
    double bth0 = 0.0;

    for(int m=0; m<ndir; ++m)
    {
    double s = 0.0;
    const int *ix = satidx.data() + size_t(m)*nw;
    const double *wx = satw.data() + size_t(m)*nw;
    for(int j=0; j<nw; ++j)
    if(ix[j]>=0)
    s += wx[j]*double(Nl[ix[j]]);
    bth[m] = s*facsat;
    bth0 = std::max(bth0,bth[m]);
    }

    const double c2 = prm.ssdsc2*sig;
    const double zcoef = c2*prm.ssdsc6, zcoefm1 = c2*(1.0-prm.ssdsc6);
    const double a0 = std::max(0.0,bth0*tmp03-prm.ssdsc4);

    for(int m=0; m<ndir; ++m)
    {
    const double a = std::max(0.0,bth[m]*tmp03-prm.ssdsc4);
    Dl[m] = -(zcoef*a0*a0 + zcoefm1*a*a);
    }
}

void seastate_st4::add(const float *N, const float *k, const float *cg, int la, int lb, int ma, int mb, double *P, double *D)
{
    const seastate_grid &g = *pg;

    for(int l=la; l<=lb; ++l)
    {
    const int b0 = g.bin(l,0);
    dissipation(N,k,cg,l,wcl.data());
        if(calmed)
        {
        for(int m=ma; m<=mb; ++m)
        D[b0+m] += wcl[m];
        continue;
        }
        for(int m=ma; m<=mb; ++m)
        {
        const int b = b0+m;
        P[b] += gam[b]*double(N[b]);
        D[b] += dst[b] + wcl[m];
        }
    }
}

// IMPHFTAIL: above the cut-off the saturation of the cut-off, F(f) = F(f_c) k_c^3 cg_c / (k^3 cg)
void seastate_st4::project(float *N, const float *k, const float *cg, int ma, int mb) const
{
    if(!prm.tail || calmed || !(ust>0.0) || lc>=nsig-1)
    return;

    const seastate_grid &g = *pg;
    const double kc = double(k[lc]);
    const double tc = kc*kc*kc*double(cg[lc]);
    const float *Nc = N + g.bin(lc,0);

    for(int l=lc+1; l<nsig; ++l)
    {
    const double kl = double(k[l]);
    const double r = tc/(kl*kl*kl*double(cg[l]))*g.sig[lc]/g.sig[l];
    float *Nl = N + g.bin(l,0);
    for(int m=ma; m<=mb; ++m)
    Nl[m] = float(r*double(Nc[m]));
    }
}
