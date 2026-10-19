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

#include"arc_sirw.h"
#include"lexer.h"
#include"ghostcell.h"
#include"slice.h"
#include<complex>
#include<cmath>
#include<iostream>

namespace
{
const double g_ = 9.81;

double kdisp(double w, double h)
{
    if(w<=0.0)
    return 0.0;
    double k = w/sqrt(g_*h);
    k = std::max(k, w*w/g_);
    for(int q=0; q<60; ++q)
    {
        const double th = tanh(k*h);
        const double dk = (g_*k*th - w*w)/(g_*th + g_*k*h*(1.0-th*th));
        k -= dk;
        if(fabs(dk)<1.0e-12*k)
        break;
    }
    return k;
}

}

void arc_sirw::ini(lexer *p, int side_, double xb, double h, double d, double dx, double flow_, double fhigh_)
{
    typedef std::complex<double> cplx;
    const cplx I(0.0,1.0);
    
    side  = side_;
    flow  = flow_;
    fhigh = fhigh_;
    hd    = h;
    x1 = xb - side*d;
    x2 = x1 - side*dx;           // x-: x2 = x1 + dx (further into the domain); x+: x2 = x1 - dx
    dtf = 1.0/(8.0*fhigh);
    
    // filter length: six periods of flow, and the travel time of the slowest band waves from the gauges to the
    // boundary with a margin (the causal response must not be cut by the taper)
    const double wh = 2.0*PI*fhigh;
    const double kh_ = kdisp(wh,h)*h;
    const double cgh = 0.5*wh/kdisp(wh,h)*(1.0 + 2.0*kh_/sinh(std::min(2.0*kh_,60.0)));
    M = int(ceil(std::max(6.0/flow, 1.5*(d+dx)/cgh)/dtf));
    const int N = 2*M+1;
    
    // frequency responses on the DFT grid f_r = r/(N dtf), r = 0..M (Hermitian for the rest): G1, G2 give the
    // wave towards the boundary at xb, the level factors omega cosh(k s_l)/sinh(kh) its velocity
    std::vector<cplx> G1(N,cplx(0.0,0.0)), G2(N,cplx(0.0,0.0));
    std::vector<double> kr(N,0.0), om(N,0.0);
    const double f0a = flow, f0b = 1.25*flow, f1a = 0.8*fhigh, f1b = fhigh;
    
    for(int r=1; r<=M; ++r)
    {
        const double f = double(r)/(double(N)*dtf);
        double w = 0.0;
        if(f>f0a && f<f1b)
        {
        w = 1.0;
        if(f<f0b) w = 0.5*(1.0-cos(PI*(f-f0a)/(f0b-f0a)));
        if(f>f1a) w = 0.5*(1.0+cos(PI*(f-f1a)/(f1b-f1a)));
        }
        if(w<=0.0)
        continue;
        
        om[r] = 2.0*PI*f;
        const double k = kdisp(om[r],h);
        kr[r] = k;
        
        // solve G1 A1 + G2 A2 = target for both wave components (rows: incident, reflected)
        // x-: target R exp(i k xb), zero for I;  x+: target I exp(-i k xb), zero for R
        const cplx a11 = exp(-I*k*x1), a12 = exp(-I*k*x2);
        const cplx a21 = exp( I*k*x1), a22 = exp( I*k*x2);
        const cplx b1v = side<0 ? cplx(0.0,0.0) : exp(-I*k*xb);
        const cplx b2v = side<0 ? exp(I*k*xb) : cplx(0.0,0.0);
        
        // det = 2 i sin(k dx): |1/(2 sin(k dx))| capped at 5 near the singularities k dx = n pi
        cplx det = a11*a22 - a12*a21;
        const double dmin = 0.2;
        if(std::abs(det)<dmin)
        det *= dmin/std::max(std::abs(det),1.0e-12);
        
        G1[r] = w*( b1v*a22 - a12*b2v)/det;
        G2[r] = w*(-b1v*a21 + a11*b2v)/det;
    }
    
    // impulse responses per level, causal part h(j dtf), j = 0..M
    // channels: 0..nlev-1 velocity per level, nlev surface elevation, nlev+1..2 nlev non-hydrostatic pressure per level
    nch = 2*nlev+1;
    h1.assign(nch,std::vector<double>(M+1,0.0));
    h2.assign(nch,std::vector<double>(M+1,0.0));
    std::vector<cplx> F1(N), F2(N);
    double ec=0.0, en=0.0;
    
    for(int l=0; l<nch; ++l)
    {
        const int lv = l>nlev ? l-nlev-1 : l;
        const double s = h*double(lv)/double(nlev-1);     // height above the bed (l = nlev: elevation channel)
        
        for(int r=0; r<N; ++r)
        F1[r] = F2[r] = cplx(0.0,0.0);
        
        for(int r=1; r<=M; ++r)
        if(kr[r]>0.0)
        {
            const double kh = std::min(kr[r]*h,30.0);
            // velocity: u = +- omega cosh(k s)/sinh(kh) eta (side: direction of travel), elevation: eta, pressure:
            // p_nh/rho = g eta (cosh(k s)/cosh(kh) - 1) (linear, p = rho g (eta - z) + p_nh)
            double fac;
            if(l<nlev)
            fac = double(side)*om[r]*(kr[r]*h<30.0 ? cosh(kr[r]*s)/sinh(kh) : exp(kr[r]*(s-h)));
            else
            if(l==nlev)
            fac = 1.0;
            else
            fac = g_*((kr[r]*h<30.0 ? cosh(kr[r]*s)/cosh(kh) : exp(kr[r]*(s-h))) - 1.0);
            F1[r] = fac*G1[r];
            F2[r] = fac*G2[r];
            F1[N-r] = std::conj(F1[r]);
            F2[N-r] = std::conj(F2[r]);
        }
        
        for(int j=-M; j<=M; ++j)
        {
            cplx s1(0.0,0.0), s2(0.0,0.0);
            for(int r=0; r<N; ++r)
            {
            const cplx ph = exp(I*(2.0*PI*double(r)*double(j)/double(N)));
            s1 += F1[r]*ph;
            s2 += F2[r]*ph;
            }
            
            if(j<0)
            {
            if(l==nlev-1)
            en += pow(s1.real(),2.0);
            continue;
            }
            
            if(l==nlev-1)
            ec += pow(s1.real(),2.0);
            
            // half-cosine taper of the last 20 % of the causal response
            double tw = 1.0;
            if(j>0.8*M)
            tw = 0.5*(1.0+cos(PI*(double(j)-0.8*M)/(0.2*M)));
            h1[l][j] = tw*s1.real()/double(N);
            h2[l][j] = tw*s2.real()/double(N);
        }
    }
    
    b1.assign(M+1,0.0);
    b2.assign(M+1,0.0);
    ycur.assign(nch,0.0);
    yprev.assign(nch,0.0);
    head = 0;
    nsample = 0;
    on = true;
    
    if(p->mpirank==0)
    std::cout<<"arc_sirw: "<<(side<0?"x- (B 98 8)":"x+ (B 99 8)")<<" gauges at x = "<<x1<<" and "<<x2<<" m, band "<<flow<<" - "<<fhigh
             <<" Hz, filters "<<M+1<<" taps at "<<dtf<<" s on "<<nlev<<" levels, non-causal energy dropped: "<<100.0*en/std::max(en+ec,1.0e-30)<<" %"<<std::endl;
}

// y-averaged surface elevation at x (linear interpolation between cell centres), all ranks
double arc_sirw::gauge(lexer *p, ghostcell *pgc, slice &eta, double xg)
{
    double sum=0.0, cnt=0.0;
    
    SLICELOOP4
    if(p->XP[IP]<=xg && xg<p->XP[IP1])
    {
        const double w = (xg-p->XP[IP])/(p->XP[IP1]-p->XP[IP]);
        sum += (1.0-w)*eta(i,j) + w*eta(i+1,j);
        cnt += 1.0;
    }
    
    sum = pgc->globalsum(sum);
    cnt = pgc->globalsum(cnt);
    
    return cnt>0.0 ? sum/cnt : 0.0;
}

void arc_sirw::update(lexer *p, ghostcell *pgc, slice &eta, double t)
{
    if(!on || p->count==count_last)
    return;
    count_last = p->count;
    
    const double e1 = gauge(p,pgc,eta,x1);
    const double e2 = gauge(p,pgc,eta,x2);
    
    if(tprev<0.0)
    {
    tprev = t;
    tnext = t;
    e1prev = e1;
    e2prev = e2;
    }
    
    // new filter samples between the last and this time step
    while(tnext<=t)
    {
        const double fac = t>tprev ? (tnext-tprev)/(t-tprev) : 1.0;
        const double s1 = e1prev + fac*(e1-e1prev);
        const double s2 = e2prev + fac*(e2-e2prev);
        
        head = (head+M)%(M+1);      // ring buffer, newest sample at head
        b1[head] = s1;
        b2[head] = s2;
        ++nsample;
        
        const int nj = std::min(nsample,M+1);
        for(int l=0; l<nch; ++l)
        {
            double y=0.0;
            for(int j=0; j<nj; ++j)
            {
            const int q = (head+j)%(M+1);
            y += h1[l][j]*b1[q] + h2[l][j]*b2[q];
            }
            yprev[l] = ycur[l];
            ycur[l]  = y;
        }
        
        tcur  = tnext;
        tnext += dtf;
    }
    
    tprev  = t;
    e1prev = e1;
    e2prev = e2;

}

double arc_sirw::velocity(double t, double z) const
{
    if(!on || tcur<0.0)
    return 0.0;
    
    const double fac = std::min(1.0,std::max(0.0,(t-tcur)/dtf));
    
    double sl = (z+hd)/hd*double(nlev-1);
    sl = std::min(double(nlev-1),std::max(0.0,sl));
    const int l0 = std::min(nlev-2,int(sl));
    const double fl = sl-double(l0);
    
    const double y0 = ycur[l0]   + fac*(ycur[l0]-yprev[l0]);
    const double y1 = ycur[l0+1] + fac*(ycur[l0+1]-yprev[l0+1]);
    
    return y0 + fl*(y1-y0);
}

double arc_sirw::pressure(double t, double z) const
{
    if(!on || tcur<0.0)
    return 0.0;
    
    const double fac = std::min(1.0,std::max(0.0,(t-tcur)/dtf));
    double sl = (z+hd)/hd*double(nlev-1);
    sl = std::min(double(nlev-1),std::max(0.0,sl));
    const int l0 = std::min(nlev-2,int(sl));
    const double fl = sl-double(l0);
    const int c0 = nlev+1+l0;
    const double y0 = ycur[c0]   + fac*(ycur[c0]-yprev[c0]);
    const double y1 = ycur[c0+1] + fac*(ycur[c0+1]-yprev[c0+1]);
    return y0 + fl*(y1-y0);
}

double arc_sirw::eta(double t) const
{
    if(!on || tcur<0.0)
    return 0.0;
    
    const double fac = std::min(1.0,std::max(0.0,(t-tcur)/dtf));
    return ycur[nlev] + fac*(ycur[nlev]-yprev[nlev]);
}

void arc_sirw::ckpt_state(std::vector<double> &v) const
{
    v.push_back(double(count_last)); v.push_back(double(head)); v.push_back(double(nsample));
    v.push_back(tnext); v.push_back(tprev); v.push_back(e1prev); v.push_back(e2prev); v.push_back(tcur);
    for(const std::vector<double> *a : {&b1,&b2,&ycur,&yprev})
    {
        v.push_back(double(a->size()));
        v.insert(v.end(),a->begin(),a->end());
    }
}

void arc_sirw::ckpt_restore(const std::vector<double> &v, size_t &k)
{
    if(k+8>v.size())
    {
        k = v.size();
        return;
    }
    count_last = int(v[k++]); head = int(v[k++]); nsample = int(v[k++]);
    tnext = v[k++]; tprev = v[k++]; e1prev = v[k++]; e2prev = v[k++]; tcur = v[k++];
    for(std::vector<double> *a : {&b1,&b2,&ycur,&yprev})
    {
        const size_t n = (k<v.size()) ? size_t(v[k++]) : 0;
        if(k+n>v.size())
        {
            k = v.size();
            return;
        }
        a->assign(v.begin()+k,v.begin()+k+n);
        k += n;
    }
}
