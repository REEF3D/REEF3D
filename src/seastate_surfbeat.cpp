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

#include"seastate_surfbeat.h"
#include"seastate_grid.h"
#include"seastate_dispersion.h"
#include<algorithm>
#include<cmath>
#include<complex>
#include<cstdint>
#include<random>

namespace
{
    const double pi = 3.14159265358979323846;

    // portable uniform [0,1) from mt19937 (identical on every platform and rank)
    double uniform(std::mt19937 &rng) {return double(rng())/4294967296.0;}
}

seastate_surfbeat::seastate_surfbeat(const seastate_grid &gin, const float *Nin, double trec_, int seed) : g(gin)
{
    const int nsig = g.nsig, ndir = g.ndir;

    // frequency spectrum E(sig), moments, T_rep = T_m-1,0 (XBeach tpDcalc)
    std::vector<double> Es(nsig,0.0);
    double mm1 = 0.0;
    m0 = 0.0;

    Dbar.assign(ndir,0.0);

    for(int l=0; l<nsig; ++l)
    {
        for(int m=0; m<ndir; ++m)
        {
        const double e = g.sig[l]*double(Nin[g.bin(l,m)]);
        Es[l] += e*g.dtheta;
        Dbar[m] += e*g.dsig[l];
        }

    m0  += Es[l]*g.dsig[l];
    mm1 += Es[l]/g.sig[l]*g.dsig[l];
    }

    if(!(m0>0.0))
    return;

    for(int m=0; m<ndir; ++m)
    Dbar[m] /= m0;

    trep = 2.0*pi*mm1/m0;
    frep = 1.0/trep;

    // frequency range: E(sig) >= 1e-3 of the peak
    const double emax = *std::max_element(Es.begin(),Es.end());
    int la = nsig-1, lb = 0;
    for(int l=0; l<nsig; ++l)
    if(Es[l]>=1.0e-3*emax)
    {
    la = std::min(la,l);
    lb = std::max(lb,l);
    }

    trec = trec_;
    df = 1.0/trec;

    const int m1 = std::max(1,int(std::ceil(g.f[la]/df)));
    const int m2 = int(std::floor(g.f[lb]/df));

    std::mt19937 rng{std::uint32_t(seed)};

    double var = 0.0;

    for(int mm=m1; mm<=m2; ++mm)
    {
    const double f = mm*df;
    const double sig = 2.0*pi*f;

    // E(sig) interpolated linearly between the grid frequencies; S(f) = 2 pi E(sig)
    int l = 0;
    while(l<nsig-2 && g.sig[l+1]<sig)
    ++l;
    const double w = std::min(std::max((sig-g.sig[l])/(g.sig[l+1]-g.sig[l]),0.0),1.0);
    const double S = 2.0*pi*((1.0-w)*Es[l] + w*Es[l+1]);

        if(!(S>0.0))
        continue;

    // direction: random draw from the directional distribution of the nearest grid frequency
    const int ln = (w<0.5) ? l : l+1;
    double sum = 0.0;
    int nz = 0;
    for(int m=0; m<ndir; ++m)
    {
    sum += double(Nin[g.bin(ln,m)]);
    if(Nin[g.bin(ln,m)]>0.0f)
    ++nz;
    }

    double th = 0.0;
    const double u = uniform(rng)*sum;
    double c = 0.0;
        for(int m=0; m<ndir; ++m)
        {
        c += double(Nin[g.bin(ln,m)]);
            if(u<c || m==ndir-1)
            {
            th = g.theta[m];
            break;
            }
        }

    const double jitter = uniform(rng);
    if(nz>1)
    th += (jitter-0.5)*g.dtheta;

    fk.push_back(f);
    ak.push_back(std::sqrt(2.0*S*df));
    thk.push_back(th);
    phk.push_back(2.0*pi*uniform(rng));

    var += 0.5*ak.back()*ak.back();
    }

    K = int(fk.size());

    // the components carry exactly the variance m0 of the spectrum
    if(var>0.0)
    for(int m=0; m<K; ++m)
    ak[m] *= std::sqrt(m0/var);

}

double seastate_surfbeat::herbers(double f1, double th1, double k1, double f2, double th2, double k2, double h, double &c3, double &th3)
{
    const double grav = seastate_gravity;
    const double w1 = 2.0*pi*f1, w2 = 2.0*pi*f2;
    const double dth = std::fabs(th2-th1) + pi;

    const double k3 = std::sqrt(std::max(k1*k1 + k2*k2 + 2.0*k1*k2*std::cos(dth),1.0e-24));
    th3 = std::atan2(k2*std::sin(th2)-k1*std::sin(th1), k2*std::cos(th2)-k1*std::cos(th1));

    c3 = 2.0*pi*(f2-f1)/k3;
    c3 = std::min(c3, 0.8*std::sqrt(grav/k3*std::tanh(k3*h)));

    const double term1 = -w1*w2;
    const double term2 = -w1 + w2;
    const double term2new = c3*k3;
    const double chk1 = std::cosh(std::min(k1*h,30.0));
    const double chk2 = std::cosh(std::min(k2*h,30.0));
    const double chk3 = std::cosh(std::min(k3*h,30.0));

    double D = -grav*k1*k2*std::cos(dth)/2.0/term1
             + grav*term2*(chk1*chk2)/((grav*k3*std::tanh(k3*h) - term2new*term2new)*term1*chk3)
             * (term2*(term1*term1/grav/grav - k1*k2*std::cos(dth))
               - 0.5*((-w1)*k2*k2/(chk2*chk2) + w2*k1*k1/(chk1*chk1)));

    // surface elevation instead of bottom pressure (Van Dongeren et al. 2003, eq. 18)
    D *= chk3/(chk1*chk2);

    return D;
}

void seastate_surfbeat::series(const std::vector<double> &y, const std::vector<double> &depth, bool bound)
{
    const int nr = int(y.size());

    nt = std::max(16,int(std::lround(trec/(0.1*trep))));
    dtbc = trec/double(nt);

    Et.assign(nr,std::vector<double>(nt,0.0));
    Zt.assign(nr,std::vector<double>(nt,0.0));
    Qx.assign(nr,std::vector<double>(nt,0.0));
    Qy.assign(nr,std::vector<double>(nt,0.0));

    if(K==0)
    return;

    const int mbase = int(std::lround(fk[0]/df));
    const int K2 = K;       // difference indices 1 .. K-1

    std::vector<std::complex<double>> W(nt);
    for(int n=0; n<nt; ++n)
    W[n] = std::polar(1.0,2.0*pi*double(n)/double(nt));

    std::vector<double> kk(K);

    for(int r=0; r<nr; ++r)
    {
    const double h = std::max(depth[r],1.0e-3);

    for(int m=0; m<K; ++m)
    kk[m] = seastate_wavenumber(2.0*pi*fk[m],h);

    // envelope: Z(t_n) = sum_m c_m exp(i 2 pi M_m n/nt), M_m the integer frequency index
    std::vector<std::complex<double>> cm(K);
    std::vector<int> im(K);
        for(int m=0; m<K; ++m)
        {
        cm[m] = std::polar(ak[m], -kk[m]*std::sin(thk[m])*y[r] + phk[m]);
        im[m] = mbase + m;
        }

        for(int n=0; n<nt; ++n)
        {
        std::complex<double> z(0.0,0.0);

            for(int m=0; m<K; ++m)
            z += cm[m]*W[(long(im[m])*n)%nt];

        Et[r][n] = 0.5*std::norm(z);
        }

        if(!bound)
        continue;

    // bound wave: difference-frequency coefficients G[d] (d = n - m)
    std::vector<std::complex<double>> G(K2,0.0), Gx(K2,0.0), Gy(K2,0.0);

        for(int m=0; m<K; ++m)
        for(int n=m+1; n<K; ++n)
        {
        const double dfr = fk[n]-fk[m];

            if(!(dfr<fk[m]))
            continue;

        double c3, th3;
        const double D = herbers(fk[m],thk[m],kk[m],fk[n],thk[n],kk[n],h,c3,th3);

        const double dky = kk[n]*std::sin(thk[n]) - kk[m]*std::sin(thk[m]);
        const std::complex<double> amp = std::polar(D*ak[m]*ak[n], -dky*y[r] + phk[n] - phk[m]);

        const int d = n-m;
        G[d]  += amp;
        Gx[d] += amp*c3*std::cos(th3);
        Gy[d] += amp*c3*std::sin(th3);
        }

        for(int n=0; n<nt; ++n)
        {
        std::complex<double> z(0.0,0.0), zx(0.0,0.0), zy(0.0,0.0);

            for(int d=1; d<K2; ++d)
            {
            const std::complex<double> e = W[(long(d)*n)%nt];
            z  += G[d]*e;
            zx += Gx[d]*e;
            zy += Gy[d]*e;
            }

        Zt[r][n] = z.real();
        Qx[r][n] = zx.real();
        Qy[r][n] = zy.real();
        }
    }
}

void seastate_surfbeat::at(int r, double t, double &E, double &zeta, double &qx, double &qy) const
{
    E = zeta = qx = qy = 0.0;

    if(nt==0 || r<0 || r>=int(Et.size()))
    return;

    double s = std::fmod(std::max(t,0.0),trec)/dtbc;
    int n0 = int(s);
    const double w = s - double(n0);
    n0 %= nt;
    const int n1 = (n0+1)%nt;

    E    = (1.0-w)*Et[r][n0] + w*Et[r][n1];
    zeta = (1.0-w)*Zt[r][n0] + w*Zt[r][n1];
    qx   = (1.0-w)*Qx[r][n0] + w*Qx[r][n1];
    qy   = (1.0-w)*Qy[r][n0] + w*Qy[r][n1];
}
