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

// FFT of the zero-padded first-order paddle signal and selection of the
// components X(t) = sum_n Re(Xc_n exp(i omega_n t)) used by the theory
int select_components(double g, double h, double fmin, double fmax, int nmax,
                      const std::vector<double> &x, double dt,
                      std::vector<cplx> &F, int &nfft, double &dw,
                      std::vector<int> &cand, std::vector<double> &Xre, std::vector<double> &Xim)
{
    const int npts = int(x.size());

    nfft=1;
    while(nfft<2*npts)
    nfft<<=1;

    F.assign(nfft,cplx(0.0,0.0));
    for(int q=0; q<npts; ++q)
    F[q] = x[q];
    fft(F,false);

    dw = 2.0*pi_/(double(nfft)*dt);

    double amax=0.0;
    for(int b=1; b<nfft/2; ++b)
    amax = std::max(amax,std::abs(F[b]));

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
        if(f>=fmin && std::abs(F[b])>=1.0e-3*amax)
        cand.push_back(b);
    }

    // keep the nmax most energetic components
    if(int(cand.size())>nmax)
    {
        std::sort(cand.begin(),cand.end(),[&](int a, int b){ return std::abs(F[a])>std::abs(F[b]); });
        cand.resize(nmax);
        std::sort(cand.begin(),cand.end());
    }

    const int N = int(cand.size());
    Xre.resize(N);
    Xim.resize(N);
    for(int n=0; n<N; ++n)
    {
        Xre[n] = 2.0*F[cand[n]].real()/double(nfft);
        Xim[n] = 2.0*F[cand[n]].imag()/double(nfft);
    }
    return N;
}

struct wm_mode
{
    int n;          // owning first-order frequency
    cplx k,ch,sh;   // wavenumber, cosh(kh), sinh(kh)
    cplx Ct;        // C*cosh(kh) = i g a/omega
    cplx U,W,P,G;   // u, w, phi_t, d/dz(phi_tt+g phi_z) at z=0
};
}

wave_lib_wavemaker2nd::wave_lib_wavemaker2nd()
{
}

wave_lib_wavemaker2nd::~wave_lib_wavemaker2nd()
{
}

void wave_lib_wavemaker2nd::compute_components(double g, double h, int shape, double zs, double ze,
                       int mode, int J, int addQ, double dw, double w2min,
                       const std::vector<int> &bin, const std::vector<double> &Xre, const std::vector<double> &Xim,
                       std::vector<double> &X2re, std::vector<double> &X2im)
{
    const int N = int(bin.size());

    // paddle shape f(s) = al + be*s on [sa,h], s = z+h
    double al=1.0, be=0.0, sa=0.0;
    if(shape==2)
    {
        al = -zs/(ze-zs);
        be = 1.0/(ze-zs);
        sa = std::max(zs,0.0);
    }

    auto If = [&](cplx k){ return Lc(k,al,be,sa,h); };   // int f cosh(k s) ds

    int bmax=0;
    for(int n=0; n<N; ++n)
    bmax = std::max(bmax,bin[n]);

    const int nb = 2*bmax+1;
    std::vector<cplx> R(nb,cplx(0.0,0.0));

    // free progressive second-order wavenumbers per bin
    std::vector<double> Kf(nb,0.0);
    for(int b=1; b<nb; ++b)
    Kf[b] = k_prog(double(b)*dw,h,g);

    // first-order modes
    std::vector<wm_mode> md;
    md.reserve(size_t(N)*size_t(J+1));
    for(int n=0; n<N; ++n)
    {
        const double w = double(bin[n])*dw;
        const cplx X(Xre[n],Xim[n]);
        for(int j=0; j<=J; ++j)
        {
            wm_mode m;
            m.n = n;
            m.k = (j==0) ? cplx(k_prog(w,h,g),0.0) : cplx(0.0,-k_evan(w,h,g,j));
            m.ch = cosh(m.k*h);
            m.sh = sinh(m.k*h);
            const cplx N0 = (2.0*m.k*h + 2.0*m.sh*m.ch)/(4.0*m.k);   // int cosh^2(k s) ds
            const cplx a  = I_*X*m.sh*If(m.k)/N0;                     // first-order amplitude
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

        // paddle terms Q = X_z phi_z - X phi_xx
        if(addQ==1)
        for(int n=0; n<N; ++n)
        {
            const cplx X(Xre[n],Xim[n]);

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

                // int f cosh(k2 s) cosh(kf s) ds and int f' sinh(k2 s) cosh(kf s) ds
                const cplx Ic = 0.5*(Lc(k2+kf,al,be,sa,h) + Lc(k2-kf,al,be,sa,h));
                const cplx Is = 0.5*be*(Ls(k2+kf,sa,h) + Ls(k2-kf,sa,h));

                cplx r = 0.5*X*(Ct2/ch2)*(-k2*k2*Ic - k2*Is);

                if(b<0)
                r = conj(r);
                R[ab] += r;
            }
        }
    }

    // paddle correction: i Omega X2 If(Kf) = R
    X2re.assign(nb,0.0);
    X2im.assign(nb,0.0);
    for(int b=1; b<nb; ++b)
    if(std::abs(R[b])>0.0)
    {
        const double Om = double(b)*dw;
        const cplx X2 = R[b]/(I_*Om*If(cplx(Kf[b],0.0)));
        X2re[b] = X2.real();
        X2im[b] = X2.imag();
    }
}

int wave_lib_wavemaker2nd::compute(double g, double h, int shape, double zs, double ze,
                       int mode, int J, int addQ, double fmin, double fmax, double f2min, int nmax,
                       const std::vector<double> &x, double dt, std::vector<double> &x2)
{
    const int npts = int(x.size());
    x2.assign(npts,0.0);
    if(npts<8)
    return 0;

    std::vector<cplx> F;
    std::vector<int> cand;
    std::vector<double> Xre, Xim, X2re, X2im;
    int nfft;
    double dw;
    const int N = select_components(g,h,fmin,fmax,nmax,x,dt,F,nfft,dw,cand,Xre,Xim);
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
        if(std::abs(F[cand[n]])>std::abs(F[bp]))
        bp=cand[n];
        f2min = 0.1*double(bp)*dw/(2.0*pi_);
    }

    compute_components(g,h,shape,zs,ze,mode,J,addQ,dw,2.0*pi_*f2min,cand,Xre,Xim,X2re,X2im);

    // synthesis on a grid with the same frequency spacing, fine enough for 2*fmax
    const int nb = int(X2re.size());
    int nsyn = nfft, r = 1;
    while(nsyn/2 <= nb)
    {
        nsyn<<=1;
        r<<=1;
    }
    std::vector<cplx> Y(nsyn,cplx(0.0,0.0));
    for(int b=1; b<nb; ++b)
    Y[b] = cplx(X2re[b],X2im[b]);
    fft(Y,true);

    const double dts = dt/double(r);
    for(int q=0; q<npts; ++q)
    {
        // sample q at t = q*dt  ->  fine index q*r
        x2[q] = Y[size_t(q)*size_t(r)].real();
    }
    (void)dts;

    return N;
}

void wave_lib_wavemaker2nd::compute_paddle_Q(double g, double h, int shape, double zs, double ze,
                       int J, double fmin, double fmax, int nmax,
                       const std::vector<double> &x, double dt, int nlev, std::vector<float> &Q)
{
    const int npts = int(x.size());
    Q.assign(size_t(nlev)*size_t(npts),0.0f);
    if(npts<8 || nlev<2)
    return;

    std::vector<cplx> F;
    std::vector<int> cand;
    std::vector<double> Xre, Xim;
    int nfft;
    double dw;
    const int N = select_components(g,h,fmin,fmax,nmax,x,dt,F,nfft,dw,cand,Xre,Xim);
    if(N==0)
    return;

    double al=1.0, be=0.0, sa=0.0;
    if(shape==2)
    {
        al = -zs/(ze-zs);
        be = 1.0/(ze-zs);
        sa = std::max(zs,0.0);
    }

    // band-limited first-order displacement X(t)
    std::vector<cplx> Y(nfft,cplx(0.0,0.0));
    for(int n=0; n<N; ++n)
    Y[cand[n]] = cplx(Xre[n],Xim[n]);
    fft(Y,true);
    std::vector<double> xs(npts);
    for(int q=0; q<npts; ++q)
    xs[q] = Y[q].real();

    // modal coefficients Ct = i g a/omega and wavenumbers
    std::vector<cplx> km, Ctm, chm;
    km.reserve(size_t(N)*(J+1));
    for(int n=0; n<N; ++n)
    {
        const double w = double(cand[n])*dw;
        const cplx X(Xre[n],Xim[n]);
        for(int j=0; j<=J; ++j)
        {
            const cplx k = (j==0) ? cplx(k_prog(w,h,g),0.0) : cplx(0.0,-k_evan(w,h,g,j));
            const cplx ch = cosh(k*h), sh = sinh(k*h);
            const cplx N0 = (2.0*k*h + 2.0*sh*ch)/(4.0*k);
            const cplx a  = I_*X*sh*Lc(k,al,be,sa,h)/N0;
            km.push_back(k);
            chm.push_back(ch);
            Ctm.push_back(I_*g*a/w);
        }
    }

    std::vector<cplx> A(nfft), B(nfft);
    for(int l=0; l<nlev; ++l)
    {
        const double s = h*double(l)/double(nlev-1);
        const double f  = (s>=sa) ? al+be*s : 0.0;
        const double fp = (s>=sa) ? be : 0.0;

        if(f==0.0 && fp==0.0)
        continue;

        std::fill(A.begin(),A.end(),cplx(0.0,0.0));
        std::fill(B.begin(),B.end(),cplx(0.0,0.0));

        for(int n=0; n<N; ++n)
        for(int j=0; j<=J; ++j)
        {
            const int m = n*(J+1)+j;
            const cplx k = km[m];
            const cplx c = Ctm[m]/chm[m];
            A[cand[n]] += -k*k*c*cosh(k*s);   // phi_xx
            B[cand[n]] +=  k*c*sinh(k*s);     // phi_z
        }
        fft(A,true);
        fft(B,true);

        for(int q=0; q<npts; ++q)
        Q[size_t(l)*npts+q] = float(fp*xs[q]*B[q].real() - f*xs[q]*A[q].real());
    }
}

int wave_lib_wavemaker2nd::resample(double **kin, int ptnum, double &t0, double &dt, std::vector<double> &xu)
{
    // uniform grid, skipping non-increasing time stamps
    std::vector<double> tt, xx;
    tt.reserve(ptnum);
    xx.reserve(ptnum);
    for(int q=0; q<ptnum; ++q)
    if(tt.empty() || kin[q][0]>tt.back())
    {
        tt.push_back(kin[q][0]);
        xx.push_back(kin[q][1]);
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
    if(ptnum<8)
    return;

    const double g = fabs(p->W22)>0.0 ? fabs(p->W22) : 9.81;
    const int addQ = (p->B119==1 && p->A10==3) ? 1 : 0;

    if(shape==1 && p->B110==1)
    {
        // restricted piston: correction computed for a full-depth piston
        zs=0.0;
        ze=h;
        if(p->mpirank==0)
        cout<<"Wave_Lib: 2nd-order correction assumes a full-depth piston (B110 ignored for the correction)"<<endl;
    }

    if(shape==2 && ze<h && p->mpirank==0)
    cout<<"Wave_Lib: WARNING flap top B111_ze below still water level, 2nd-order correction not reliable"<<endl;

    double t0, dt;
    std::vector<double> xu, x2;
    const int nu = resample(kin,ptnum,t0,dt,xu);
    if(nu<8)
    return;

    const double starttime = pgc->timer();

    const int N = compute(g,h,shape,zs,ze,p->B113,p->B113_J,addQ,p->B114_fmin,p->B114_fmax,p->B114_f2min,500,xu,dt,x2);

    // interpolate the correction back onto the input time stamps
    double x2max=0.0, xmax=0.0;
    for(int q=0; q<ptnum; ++q)
    {
        const double s = (kin[q][0]-t0)/dt;
        int iq = int(floor(s));
        double corr=0.0;
        if(iq>=0 && iq<nu-1)
        corr = x2[iq] + (s-double(iq))*(x2[iq+1]-x2[iq]);
        else
        if(iq>=nu-1)
        corr = x2[nu-1];

        xmax  = std::max(xmax,fabs(kin[q][1]));
        x2max = std::max(x2max,fabs(corr));
        kin[q][1] += corr;
    }

    if(p->mpirank==0)
    {
        cout<<"Wave_Lib: 2nd-order wavemaker correction ";
        if(p->B113==1) cout<<"(sub+super)";
        if(p->B113==2) cout<<"(subharmonic)";
        if(p->B113==3) cout<<"(superharmonic)";
        cout<<"  components: "<<N<<"  evanescent modes: "<<p->B113_J<<"  paddle BC terms: "<<addQ<<endl;
        cout<<"Wave_Lib: max |X1|: "<<xmax<<"  max |X2|: "<<x2max<<"  time: "<<pgc->timer()-starttime<<" s"<<endl;
    }
}

void wave_lib_wavemaker2nd::make_Qtable(lexer *p, ghostcell *pgc, double **kin, int ptnum, double h, int shape, double zs, double ze)
{
    const double g = fabs(p->W22)>0.0 ? fabs(p->W22) : 9.81;

    if(shape==1)
    {
        zs=0.0;
        ze=h;
    }

    std::vector<double> xu;
    const int nu = resample(kin,ptnum,Qt0,Qdt,xu);
    if(nu<8)
    return;

    const double starttime = pgc->timer();

    Qh = h;
    Qnlev = 41;
    Qnt = nu;
    compute_paddle_Q(g,h,shape,zs,ze,p->B113_J,p->B114_fmin,p->B114_fmax,500,xu,Qdt,Qnlev,Qtab);

    if(p->mpirank==0)
    cout<<"Wave_Lib: moving-paddle BC terms tabulated, levels: "<<Qnlev<<"  samples: "<<Qnt<<"  evanescent modes: "<<p->B113_J<<"  time: "<<pgc->timer()-starttime<<" s"<<endl;
}
#endif
