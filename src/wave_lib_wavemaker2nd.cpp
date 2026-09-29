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

#include"wave_lib_wavemaker2nd.h"
#ifndef WM2ND_STANDALONE
#include"lexer.h"
#include"ghostcell.h"
#endif
#include<complex>
#include<cmath>
#include<algorithm>
#include<iostream>

// Second-order wavemaker theory following Schaffer (1996, Ocean Eng. 23(1)),
// in the notation of Akselsen (2025, arXiv:2502.09586), eqs. (3), (15), (19)-(23).
// z-integrals are evaluated in closed form, the first-order solution contains
// the progressive mode and J evanescent modes per frequency.

namespace
{
typedef std::complex<double> cplx;
const double pi_ = 3.14159265358979323846;
const cplx I_(0.0,1.0);

void fft(std::vector<cplx> &a, bool inverse)   // unnormalized; forward exp(-i..), inverse exp(+i..)
{
    const int n = int(a.size());
    for(int i=1, j=0; i<n; ++i)
    {
        int bit = n>>1;
        for(; j&bit; bit>>=1)
        j ^= bit;
        j ^= bit;
        if(i<j)
        std::swap(a[i],a[j]);
    }
    for(int len=2; len<=n; len<<=1)
    {
        const double ang = 2.0*pi_/double(len)*(inverse?1.0:-1.0);
        for(int j=0; j<len/2; ++j)
        {
            const cplx w = std::polar(1.0,ang*double(j));
            for(int i=0; i<n; i+=len)
            {
                cplx u = a[i+j];
                cplx v = a[i+j+len/2]*w;
                a[i+j] = u+v;
                a[i+j+len/2] = u-v;
            }
        }
    }
}

double k_prog(double w, double h, double g)
{
    if(w<=0.0)
    return 0.0;
    // Fenton & McKee initial guess, then Newton
    const double alpha = w*w*h/g;
    double k = alpha*pow(tanh(pow(alpha,0.75)),-2.0/3.0)/h;
    for(int it=0; it<100; ++it)
    {
        const double th = tanh(k*h);
        const double F  = g*k*th - w*w;
        const double dF = g*th + g*k*h*(1.0-th*th);
        const double dk = F/dF;
        k -= dk;
        if(fabs(dk)<1.0e-15*k)
        break;
    }
    return k;
}

// j-th evanescent root: w^2 = -g kap tan(kap h), kap h in ((j-1/2)pi, j pi)
double k_evan(double w, double h, double g, int j)
{
    double lo = (double(j)-0.5)*pi_/h*(1.0+1.0e-14);
    double hi = double(j)*pi_/h;
    for(int it=0; it<200; ++it)
    {
        const double mid = 0.5*(lo+hi);
        const double F = w*w + g*mid*tan(mid*h);
        if(F>0.0)
        hi=mid;
        else
        lo=mid;
        if(hi-lo<1.0e-15*hi)
        break;
    }
    return 0.5*(lo+hi);
}

// int_a^b sinh(q s) ds
cplx Ls(cplx q, double a, double b)
{
    if(std::abs(q)*b<1.0e-7)
    return q*(b*b-a*a)*0.5;

    return 2.0*sinh(0.5*q*(b+a))*sinh(0.5*q*(b-a))/q;
}

// int_a^b (al + be s) cosh(q s) ds
cplx Lc(cplx q, double al, double be, double a, double b)
{
    if(std::abs(q)*b<1.0e-7)
    return al*(b-a) + 0.5*be*(b*b-a*a);

    return ((al+be*b)*sinh(q*b) - (al+be*a)*sinh(q*a))/q - be*Ls(q,a,b)/q;
}

// sinh(x h)/x
cplx sinhc(cplx x, double h)
{
    if(std::abs(x)*h<1.0e-7)
    return cplx(h,0.0);

    return sinh(x*h)/x;
}

struct wm_mode
{
    int n;          // owning first-order frequency
    cplx k,ch,sh;   // wavenumber, cosh(kh), sinh(kh)
    cplx Ct;        // C*cosh(kh) = i g a/omega
    cplx U,W,P,G;   // u, w, phi_t, d/dz(phi_tt+g phi_z) at z=0
};

typedef wave_lib_wavemaker2nd::shape_t shape_t;

// int f_p cosh(k s) ds over [sa,h]
cplx If_shape(const shape_t &s, cplx k, double h)
{
    return Lc(k,s.al,s.be,s.sa,h);
}

// FFT of the zero-padded first-order paddle signals and selection of the
// components X_p(t) = sum_n Re(X[p][n] exp(i omega_n t)) used by the theory
int select_components(double g, double h, double fmin, double fmax, int nmax,
                      const std::vector<std::vector<double> > &x, double dt,
                      int &nfft, double &dw, std::vector<int> &cand,
                      std::vector<std::vector<cplx> > &X, std::vector<double> &amp)
{
    const int P = int(x.size());
    const int npts = int(x[0].size());

    nfft=1;
    while(nfft<2*npts)
    nfft<<=1;

    std::vector<std::vector<cplx> > F(P);
    for(int p=0; p<P; ++p)
    {
        F[p].assign(nfft,cplx(0.0,0.0));
        for(int q=0; q<npts; ++q)
        F[p][q] = x[p][q];
        fft(F[p],false);
    }

    dw = 2.0*pi_/(double(nfft)*dt);

    // combined amplitude of all segments
    amp.assign(nfft/2,0.0);
    double amax=0.0;
    for(int b=1; b<nfft/2; ++b)
    {
        double a2=0.0;
        for(int p=0; p<P; ++p)
        a2 += std::norm(F[p][b]);
        amp[b] = sqrt(a2);
        amax = std::max(amax,amp[b]);
    }

    cand.clear();
    for(int b=1; b<nfft/2; ++b)
    {
        const double f = double(b)*dw/(2.0*pi_);
        const double kh = k_prog(double(b)*dw,h,g)*h;
        if(kh>100.0)
        continue;

        if(fmax>0.0)
        {
            if(f>=fmin && f<=fmax)
            cand.push_back(b);
        }
        else
        if(f>=fmin && amp[b]>=1.0e-3*amax)
        cand.push_back(b);
    }

    // keep the nmax most energetic components
    if(int(cand.size())>nmax)
    {
        std::sort(cand.begin(),cand.end(),[&](int a, int b){ return amp[a]>amp[b]; });
        cand.resize(nmax);
        std::sort(cand.begin(),cand.end());
    }

    const int N = int(cand.size());
    X.assign(P,std::vector<cplx>(N));
    for(int p=0; p<P; ++p)
    for(int n=0; n<N; ++n)
    X[p][n] = 2.0*F[p][cand[n]]/double(nfft);

    return N;
}

shape_t make_shape(int shape, double zs, double ze)
{
    if(shape==2)
    return wave_lib_wavemaker2nd::shape_flap(zs,ze);

    return wave_lib_wavemaker2nd::shape_piston();
}
}

wave_lib_wavemaker2nd::shape_t wave_lib_wavemaker2nd::shape_piston()
{
    shape_t s;
    s.al = 1.0;
    s.be = 0.0;
    s.sa = 0.0;
    return s;
}

wave_lib_wavemaker2nd::shape_t wave_lib_wavemaker2nd::shape_flap(double zs, double ze)
{
    shape_t s;
    s.al = -zs/(ze-zs);
    s.be = 1.0/(ze-zs);
    s.sa = std::max(zs,0.0);
    return s;
}

wave_lib_wavemaker2nd::wave_lib_wavemaker2nd()
{
}

wave_lib_wavemaker2nd::~wave_lib_wavemaker2nd()
{
}

void wave_lib_wavemaker2nd::compute_components_multi(double g, double h, const std::vector<shape_t> &shp,
                       int mode, int J, int addQ, double dw, double w2min,
                       const std::vector<int> &bin, const std::vector<std::vector<cplx> > &X,
                       std::vector<std::vector<cplx> > &X2, std::vector<cplx> *RQ)
{
    const int N = int(bin.size());
    const int P = int(shp.size());

    int bmax=0;
    for(int n=0; n<N; ++n)
    bmax = std::max(bmax,bin[n]);

    const int nb = 2*bmax+1;
    std::vector<cplx> R(nb,cplx(0.0,0.0));

    // free progressive second-order wavenumbers per bin
    std::vector<double> Kf(nb,0.0);
    for(int b=1; b<nb; ++b)
    Kf[b] = k_prog(double(b)*dw,h,g);

    // first-order modes of the combined paddle
    std::vector<wm_mode> md;
    md.reserve(size_t(N)*size_t(J+1));
    for(int n=0; n<N; ++n)
    {
        const double w = double(bin[n])*dw;
        for(int j=0; j<=J; ++j)
        {
            wm_mode m;
            m.n = n;
            m.k = (j==0) ? cplx(k_prog(w,h,g),0.0) : cplx(0.0,-k_evan(w,h,g,j));
            m.ch = cosh(m.k*h);
            m.sh = sinh(m.k*h);
            const cplx N0 = (2.0*m.k*h + 2.0*m.sh*m.ch)/(4.0*m.k);   // int cosh^2(k s) ds
            cplx XI(0.0,0.0);
            for(int p=0; p<P; ++p)
            XI += X[p][n]*If_shape(shp[p],m.k,h);
            const cplx a  = I_*m.sh*XI/N0;                            // first-order amplitude
            const cplx th = w*w/(g*m.k);                              // tanh(kh) from dispersion
            m.Ct = I_*g*a/w;
            m.U = -I_*m.k*m.Ct;
            m.W = m.k*m.Ct*th;
            m.P = I_*w*m.Ct;
            m.G = m.Ct*m.k*(g*m.k - w*w*th);
            md.push_back(m);
        }
    }
    const int M = int(md.size());

    for(int sgn=1; sgn>=-1; sgn-=2)
    {
        if(sgn==1 && mode==2)
        continue;
        if(sgn==-1 && mode==3)
        continue;

        // bound wave contribution (-phi^(21)_x at the paddle, projected)
        if(RQ==nullptr)
        for(int q1=0; q1<M; ++q1)
        {
            const wm_mode &a1 = md[q1];
            const double w1 = double(bin[a1.n])*dw;

            for(int q2=0; q2<M; ++q2)
            {
                const wm_mode &a2 = md[q2];
                int b = bin[a1.n] + sgn*bin[a2.n];
                if(double(std::abs(b))*dw<w2min || b==0)
                continue;

                const double Om = w1 + double(sgn)*double(bin[a2.n])*dw;

                cplx k2=a2.k, ch2=a2.ch, sh2=a2.sh, U2=a2.U, W2=a2.W, G2=a2.G;
                if(sgn==-1)
                {
                    k2=conj(k2); ch2=conj(ch2); sh2=conj(sh2); U2=conj(U2); W2=conj(W2); G2=conj(G2);
                }

                const cplx K   = a1.k + k2;
                const cplx chK = a1.ch*ch2 + a1.sh*sh2;
                const cplx shK = a1.sh*ch2 + a1.ch*sh2;

                const cplx D = -I_*Om*(a1.U*U2 + a1.W*W2) + a1.P*G2/g;
                const cplx B = D/(g*K*shK - Om*Om*chK);

                const int ab = std::abs(b);
                const double kf = Kf[ab];
                const double chf = cosh(kf*h), shf = sinh(kf*h);

                // int_0^h cosh(K s) cosh(kf s) ds
                const cplx Sp = (shK*chf + chK*shf);           // sinh((K+kf)h)
                cplx P1 = (std::abs(K+kf)*h<1.0e-7) ? cplx(h,0.0) : Sp/(K+kf);
                cplx P2;
                if(std::abs(K-kf)*h<0.1)
                P2 = sinhc(K-kf,h);
                else
                P2 = (shK*chf - chK*shf)/(K-kf);

                cplx r = 0.5*(-I_*K)*B*0.5*(P1+P2);

                if(b<0)
                r = conj(r);
                R[ab] += r;
            }
        }

        // paddle terms Q = X_z phi_z - X phi_xx, S = sum_p X_p f_p
        if(addQ==1)
        for(int n=0; n<N; ++n)
        {
            for(int q2=0; q2<M; ++q2)
            {
                const wm_mode &a2 = md[q2];
                int b = bin[n] + sgn*bin[a2.n];
                if(double(std::abs(b))*dw<w2min || b==0)
                continue;

                cplx k2=a2.k, ch2=a2.ch, Ct2=a2.Ct;
                if(sgn==-1)
                {
                    k2=conj(k2); ch2=conj(ch2); Ct2=conj(Ct2);
                }

                const int ab = std::abs(b);
                const double kf = Kf[ab];

                cplx sum(0.0,0.0);
                for(int p=0; p<P; ++p)
                {
                    const shape_t &s = shp[p];
                    // int f cosh(k2 s) cosh(kf s) ds and int f' sinh(k2 s) cosh(kf s) ds
                    const cplx Ic = 0.5*(Lc(k2+kf,s.al,s.be,s.sa,h) + Lc(k2-kf,s.al,s.be,s.sa,h));
                    const cplx Is = 0.5*s.be*(Ls(k2+kf,s.sa,h) + Ls(k2-kf,s.sa,h));
                    sum += X[p][n]*(-k2*k2*Ic - k2*Is);
                }

                cplx r = 0.5*(Ct2/ch2)*sum;

                if(b<0)
                r = conj(r);
                R[ab] += r;
            }
        }
    }

    // only the projected paddle terms requested (B119 forcing)
    if(RQ!=nullptr)
    {
        *RQ = R;
        X2.assign(P,std::vector<cplx>(nb,cplx(0.0,0.0)));
        return;
    }

    // direction of the paddle motion per first-order component, e = X/|X|
    // with the phase of the dominant segment removed (P=1: e=1)
    std::vector<std::vector<cplx> > e(N,std::vector<cplx>(P,cplx(1.0,0.0)));
    std::vector<int> rel;
    if(P>1)
    {
        double wmax=0.0;
        std::vector<double> wn(N,0.0);
        for(int n=0; n<N; ++n)
        {
            int pd=0;
            for(int p=0; p<P; ++p)
            {
                wn[n] += std::norm(X[p][n]);
                if(std::abs(X[p][n])>std::abs(X[pd][n]))
                pd=p;
            }
            wmax = std::max(wmax,wn[n]);
            const cplx ph = (std::abs(X[pd][n])>0.0) ? conj(X[pd][n])/std::abs(X[pd][n]) : cplx(1.0,0.0);
            const double nrm = sqrt(wn[n]);
            for(int p=0; p<P; ++p)
            e[n][p] = (nrm>0.0) ? X[p][n]*ph/nrm : cplx(p==0?1.0:0.0,0.0);
        }
        // reliable components for the interpolation: >= 5% of the peak energy
        for(int n=0; n<N; ++n)
        if(wn[n]>=0.05*wmax)
        rel.push_back(n);
    }

    // paddle correction: i Omega sum_p X2_p If_p(Kf) = R, X2_p = c e_p(Omega)
    X2.assign(P,std::vector<cplx>(nb,cplx(0.0,0.0)));
    for(int b=1; b<nb; ++b)
    if(std::abs(R[b])>0.0)
    {
        const double Om = double(b)*dw;

        std::vector<cplx> eb(P,cplx(1.0,0.0));
        if(P>1)
        {
            // linear interpolation of the direction in frequency, clamped at the band edges
            const int nr = int(rel.size());
            if(double(b)<=double(bin[rel[0]]))
            eb = e[rel[0]];
            else
            if(double(b)>=double(bin[rel[nr-1]]))
            eb = e[rel[nr-1]];
            else
            {
                int r1=1;
                while(bin[rel[r1]]<b)
                ++r1;
                const int n0=rel[r1-1], n1=rel[r1];
                const double fac = double(b-bin[n0])/double(bin[n1]-bin[n0]);
                double nrm=0.0;
                for(int p=0; p<P; ++p)
                {
                    eb[p] = e[n0][p] + fac*(e[n1][p]-e[n0][p]);
                    nrm += std::norm(eb[p]);
                }
                nrm = sqrt(nrm);
                for(int p=0; p<P; ++p)
                eb[p] /= (nrm>0.0?nrm:1.0);
            }
        }

        cplx den(0.0,0.0);
        for(int p=0; p<P; ++p)
        den += eb[p]*If_shape(shp[p],cplx(Kf[b],0.0),h);
        den *= I_*Om;

        if(std::abs(den)<1.0e-14*h)
        continue;

        const cplx c = R[b]/den;
        for(int p=0; p<P; ++p)
        X2[p][b] = c*eb[p];
    }
}

int wave_lib_wavemaker2nd::compute_multi(double g, double h, const std::vector<shape_t> &shp,
                       int mode, int J, int addQ, double fmin, double fmax, double f2min, int nmax,
                       const std::vector<std::vector<double> > &x, double dt, std::vector<std::vector<double> > &x2)
{
    const int P = int(shp.size());
    const int npts = int(x[0].size());
    x2.assign(P,std::vector<double>(npts,0.0));
    if(npts<8)
    return 0;

    std::vector<int> cand;
    std::vector<std::vector<cplx> > X, X2;
    std::vector<double> amp;
    int nfft;
    double dw;
    const int N = select_components(g,h,fmin,fmax,nmax,x,dt,nfft,dw,cand,X,amp);
    if(N==0)
    return 0;

    // lowest second-order frequency that is corrected: very long difference
    // frequencies (set-down under the ramped wave train, spectral leakage)
    // would require an ever growing paddle drift and are not corrected.
    // default: 0.1 x peak frequency of the first-order paddle signal
    if(f2min<=0.0)
    {
        int bp=cand[0];
        for(int n=0; n<N; ++n)
        if(amp[cand[n]]>amp[bp])
        bp=cand[n];
        f2min = 0.1*double(bp)*dw/(2.0*pi_);
    }

    compute_components_multi(g,h,shp,mode,J,addQ,dw,2.0*pi_*f2min,cand,X,X2);

    // synthesis on a grid with the same frequency spacing, fine enough for 2*fmax
    const int nb = int(X2[0].size());
    int nsyn = nfft, r = 1;
    while(nsyn/2 <= nb)
    {
        nsyn<<=1;
        r<<=1;
    }
    for(int p=0; p<P; ++p)
    {
        std::vector<cplx> Y(nsyn,cplx(0.0,0.0));
        for(int b=1; b<nb; ++b)
        Y[b] = X2[p][b];
        fft(Y,true);

        for(int q=0; q<npts; ++q)
        x2[p][q] = Y[size_t(q)*size_t(r)].real();
    }

    return N;
}

void wave_lib_wavemaker2nd::compute_paddle_Q_multi(double g, double h, const std::vector<shape_t> &shp,
                       int J, double fmin, double fmax, int nmax,
                       const std::vector<std::vector<double> > &x, double dt, int nlev, std::vector<float> &Q)
{
    // The pointwise Taylor terms Q(z,t) of the linear paddle solution are
    // singular at the paddle/free-surface corner (phi_xx ~ log r) and, applied
    // node by node, can exceed the first-order paddle velocity near the surface.
    // Only their projection onto the free progressive mode generates a
    // propagating second-order wave (this is also what B113 cancels), so Q is
    // replaced by that projection for every second-order frequency:
    //   Q_prog(s,t) = sum_b Re( Qb cosh(K_b s) exp(i Omega_b t) ),
    //   Qb = int Q cosh(K_b s) ds / int cosh^2(K_b s) ds.
    // It is smooth, bounded and radiates exactly the same free wave.

    const int npts = int(x[0].size());
    Q.assign(size_t(nlev)*size_t(npts),0.0f);
    if(npts<8 || nlev<2)
    return;

    std::vector<int> cand;
    std::vector<std::vector<cplx> > X, X2;
    std::vector<double> amp;
    std::vector<cplx> RQ;
    int nfft;
    double dw;
    const int N = select_components(g,h,fmin,fmax,nmax,x,dt,nfft,dw,cand,X,amp);
    if(N==0)
    return;

    compute_components_multi(g,h,shp,1,J,1,dw,0.0,cand,X,X2,&RQ);

    const int nb = int(RQ.size());
    std::vector<double> Kb(nb,0.0);
    std::vector<cplx> Qb(nb,cplx(0.0,0.0));
    for(int b=1; b<nb; ++b)
    if(std::abs(RQ[b])>0.0)
    {
        Kb[b] = k_prog(double(b)*dw,h,g);
        Qb[b] = -RQ[b];                                             // R accumulates -int Q cosh
    }

    // cosh(K s) / int_0^h cosh^2(K s) ds, overflow-safe for large K h
    auto shapefac = [&](double K, double s)
    {
        if(K*h<20.0)
        return cosh(K*s)/(0.5*h + sinh(2.0*K*h)/(4.0*K));

        return 4.0*K*exp(K*(s-2.0*h))*0.5*(1.0+exp(-2.0*K*s))/(1.0 + 2.0*K*h*exp(-2.0*K*h));
    };

    // synthesis grid with the same frequency spacing, fine enough for 2*fmax
    int nsyn = nfft, r = 1;
    while(nsyn/2 <= nb)
    {
        nsyn<<=1;
        r<<=1;
    }

    std::vector<cplx> Y(nsyn);
    for(int l=0; l<nlev; ++l)
    {
        const double s = h*double(l)/double(nlev-1);

        std::fill(Y.begin(),Y.end(),cplx(0.0,0.0));
        for(int b=1; b<nb; ++b)
        if(std::abs(Qb[b])>0.0)
        Y[b] = Qb[b]*shapefac(Kb[b],s);
        fft(Y,true);

        for(int q=0; q<npts; ++q)
        Q[size_t(l)*npts+q] = float(Y[size_t(q)*size_t(r)].real());
    }
}

// ---- single-segment wrappers ----

int wave_lib_wavemaker2nd::compute(double g, double h, int shape, double zs, double ze,
                       int mode, int J, int addQ, double fmin, double fmax, double f2min, int nmax,
                       const std::vector<double> &x, double dt, std::vector<double> &x2)
{
    std::vector<shape_t> shp(1,make_shape(shape,zs,ze));
    std::vector<std::vector<double> > xv(1,x), x2v;
    const int N = compute_multi(g,h,shp,mode,J,addQ,fmin,fmax,f2min,nmax,xv,dt,x2v);
    x2 = x2v[0];
    return N;
}

void wave_lib_wavemaker2nd::compute_components(double g, double h, int shape, double zs, double ze,
                       int mode, int J, int addQ, double dw, double w2min,
                       const std::vector<int> &bin, const std::vector<double> &Xre, const std::vector<double> &Xim,
                       std::vector<double> &X2re, std::vector<double> &X2im)
{
    std::vector<shape_t> shp(1,make_shape(shape,zs,ze));
    std::vector<std::vector<cplx> > X(1,std::vector<cplx>(bin.size())), X2;
    for(size_t n=0; n<bin.size(); ++n)
    X[0][n] = cplx(Xre[n],Xim[n]);

    compute_components_multi(g,h,shp,mode,J,addQ,dw,w2min,bin,X,X2);

    X2re.resize(X2[0].size());
    X2im.resize(X2[0].size());
    for(size_t b=0; b<X2[0].size(); ++b)
    {
        X2re[b] = X2[0][b].real();
        X2im[b] = X2[0][b].imag();
    }
}

void wave_lib_wavemaker2nd::compute_paddle_Q(double g, double h, int shape, double zs, double ze,
                       int J, double fmin, double fmax, int nmax,
                       const std::vector<double> &x, double dt, int nlev, std::vector<float> &Q)
{
    std::vector<shape_t> shp(1,make_shape(shape,zs,ze));
    std::vector<std::vector<double> > xv(1,x);
    compute_paddle_Q_multi(g,h,shp,J,fmin,fmax,nmax,xv,dt,nlev,Q);
}

int wave_lib_wavemaker2nd::resample(double **kin, int ptnum, int col, double &t0, double &dt, std::vector<double> &xu)
{
    // uniform grid, skipping non-increasing time stamps
    std::vector<double> tt, xx;
    tt.reserve(ptnum);
    xx.reserve(ptnum);
    for(int q=0; q<ptnum; ++q)
    if(tt.empty() || kin[q][0]>tt.back())
    {
        tt.push_back(kin[q][0]);
        xx.push_back(kin[q][col]);
    }

    const int nt = int(tt.size());
    if(nt<8)
    return 0;

    std::vector<double> dts(nt-1);
    for(int q=0; q<nt-1; ++q)
    dts[q] = tt[q+1]-tt[q];
    std::nth_element(dts.begin(),dts.begin()+dts.size()/2,dts.end());
    dt = dts[dts.size()/2];
    t0 = tt.front();

    const int nu = int((tt.back()-tt.front())/dt) + 1;
    xu.resize(nu);
    int qq=0;
    for(int q=0; q<nu; ++q)
    {
        const double t = t0 + double(q)*dt;
        while(qq<nt-2 && tt[qq+1]<t)
        ++qq;
        const double fac = std::min(1.0,std::max(0.0,(t-tt[qq])/(tt[qq+1]-tt[qq])));
        xu[q] = xx[qq] + fac*(xx[qq+1]-xx[qq]);
    }
    return nu;
}

double wave_lib_wavemaker2nd::paddle_Q(double t, double z) const
{
    if(Qnt<2 || Qnlev<2)
    return 0.0;

    // levels s = z+h in [0,h]; above still water the top level is used
    double sl = (z+Qh)/Qh*double(Qnlev-1);
    sl = std::min(double(Qnlev-1),std::max(0.0,sl));
    int l0 = std::min(Qnlev-2,int(sl));
    double fl = sl-double(l0);

    double st = (t-Qt0)/Qdt;
    if(st<0.0 || st>double(Qnt-1))
    return 0.0;
    int q0 = std::min(Qnt-2,int(st));
    double ft = st-double(q0);

    const float *a = &Qtab[size_t(l0)*Qnt];
    const float *b = &Qtab[size_t(l0+1)*Qnt];
    double qa = a[q0] + ft*(a[q0+1]-a[q0]);
    double qb = b[q0] + ft*(b[q0+1]-b[q0]);

    return qa + fl*(qb-qa);
}

#ifndef WM2ND_STANDALONE
void wave_lib_wavemaker2nd::correct(lexer *p, ghostcell *pgc, double **kin, int ptnum, double h, int shape, double zs, double ze)
{
    if(shape==1 && p->B110==1)
    {
        // restricted piston: correction computed for a full-depth piston
        if(p->mpirank==0)
        cout<<"Wave_Lib: 2nd-order correction assumes a full-depth piston (B110 ignored for the correction)"<<endl;
    }

    if(shape==2 && ze<h && p->mpirank==0)
    cout<<"Wave_Lib: WARNING flap top B111_ze below still water level, 2nd-order correction not reliable"<<endl;

    std::vector<shape_t> shp(1,make_shape(shape,zs,ze));
    std::vector<int> col(1,1);
    correct_impl(p,pgc,kin,ptnum,h,shp,col);
}

void wave_lib_wavemaker2nd::correct_double(lexer *p, ghostcell *pgc, double **kin, int ptnum, double h, double zs, double z2, double ze)
{
    if(ze<h && p->mpirank==0)
    cout<<"Wave_Lib: WARNING double flap top B112_ze below still water level, 2nd-order correction not reliable"<<endl;

    std::vector<shape_t> shp;
    shp.push_back(shape_flap(zs,z2));
    shp.push_back(shape_flap(z2,ze));
    std::vector<int> col;
    col.push_back(1);
    col.push_back(2);
    correct_impl(p,pgc,kin,ptnum,h,shp,col);
}

void wave_lib_wavemaker2nd::correct_impl(lexer *p, ghostcell *pgc, double **kin, int ptnum, double h,
                      const std::vector<shape_t> &shp, const std::vector<int> &col)
{
    if(ptnum<8)
    return;

    const double g = fabs(p->W22)>0.0 ? fabs(p->W22) : 9.81;
    const int addQ = (p->B119==1 && p->A10==3) ? 1 : 0;
    const int P = int(shp.size());

    double t0=0.0, dt=1.0;
    std::vector<std::vector<double> > xu(P), x2;
    int nu=0;
    for(int c=0; c<P; ++c)
    nu = resample(kin,ptnum,col[c],t0,dt,xu[c]);
    if(nu<8)
    return;

    const double starttime = pgc->timer();

    const int N = compute_multi(g,h,shp,p->B113,p->B113_J,addQ,p->B114_fmin,p->B114_fmax,p->B114_f2min,500,xu,dt,x2);

    // interpolate the correction back onto the input time stamps
    std::vector<double> x2max(P,0.0), xmax(P,0.0);
    for(int q=0; q<ptnum; ++q)
    {
        const double s = (kin[q][0]-t0)/dt;
        int iq = int(floor(s));
        for(int c=0; c<P; ++c)
        {
            double corr=0.0;
            if(iq>=0 && iq<nu-1)
            corr = x2[c][iq] + (s-double(iq))*(x2[c][iq+1]-x2[c][iq]);
            else
            if(iq>=nu-1)
            corr = x2[c][nu-1];

            xmax[c]  = std::max(xmax[c],fabs(kin[q][col[c]]));
            x2max[c] = std::max(x2max[c],fabs(corr));
            kin[q][col[c]] += corr;
        }
    }

    if(p->mpirank==0)
    {
        cout<<"Wave_Lib: 2nd-order wavemaker correction ";
        if(p->B113==1) cout<<"(sub+super)";
        if(p->B113==2) cout<<"(subharmonic)";
        if(p->B113==3) cout<<"(superharmonic)";
        cout<<"  components: "<<N<<"  evanescent modes: "<<p->B113_J<<"  paddle BC terms: "<<addQ<<endl;
        for(int c=0; c<P; ++c)
        cout<<"Wave_Lib: "<<(P>1?(c==0?"lower flap ":"upper flap "):"")<<"max |X1|: "<<xmax[c]<<"  max |X2|: "<<x2max[c]<<endl;
        cout<<"Wave_Lib: time: "<<pgc->timer()-starttime<<" s"<<endl;
    }
}

void wave_lib_wavemaker2nd::make_Qtable(lexer *p, ghostcell *pgc, double **kin, int ptnum, double h, int shape, double zs, double ze)
{
    std::vector<shape_t> shp(1,make_shape(shape,zs,ze));
    std::vector<int> col(1,1);
    make_Qtable_impl(p,pgc,kin,ptnum,h,shp,col);
}

void wave_lib_wavemaker2nd::make_Qtable_double(lexer *p, ghostcell *pgc, double **kin, int ptnum, double h, double zs, double z2, double ze)
{
    std::vector<shape_t> shp;
    shp.push_back(shape_flap(zs,z2));
    shp.push_back(shape_flap(z2,ze));
    std::vector<int> col;
    col.push_back(1);
    col.push_back(2);
    make_Qtable_impl(p,pgc,kin,ptnum,h,shp,col);
}

void wave_lib_wavemaker2nd::make_Qtable_impl(lexer *p, ghostcell *pgc, double **kin, int ptnum, double h,
                      const std::vector<shape_t> &shp, const std::vector<int> &col)
{
    const double g = fabs(p->W22)>0.0 ? fabs(p->W22) : 9.81;
    const int P = int(shp.size());

    std::vector<std::vector<double> > xu(P);
    int nu=0;
    for(int c=0; c<P; ++c)
    nu = resample(kin,ptnum,col[c],Qt0,Qdt,xu[c]);
    if(nu<8)
    return;

    const double starttime = pgc->timer();

    Qh = h;
    Qnlev = 41;
    Qnt = nu;
    compute_paddle_Q_multi(g,h,shp,p->B113_J,p->B114_fmin,p->B114_fmax,500,xu,Qdt,Qnlev,Qtab);

    if(p->mpirank==0)
    cout<<"Wave_Lib: moving-paddle BC terms tabulated, levels: "<<Qnlev<<"  samples: "<<Qnt<<"  evanescent modes: "<<p->B113_J<<"  time: "<<pgc->timer()-starttime<<" s"<<endl;
}
#endif
