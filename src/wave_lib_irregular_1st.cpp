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

#include"wave_lib_irregular_1st.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<cmath>

// Expose glibc's libmvec vector variants (normally only with -ffast-math).
// Used only when compiled with -fopenmp-simd; otherwise these are no-ops.
#if defined(_OPENMP) || defined(REEF3D_SIMD_MATH)
extern "C" {
#pragma omp declare simd notinbranch
double sin(double) noexcept;
#pragma omp declare simd notinbranch
double cos(double) noexcept;
#pragma omp declare simd notinbranch
double cosh(double) noexcept;
#pragma omp declare simd notinbranch
double exp(double) noexcept;
}
#endif

wave_lib_irregular_1st::wave_lib_irregular_1st(lexer *p, ghostcell *pgc) : wave_lib_parameters(p,pgc) 
{ 
    parameters(p,pgc);
    
	if(p->B85!=4 && p->B85!=5 && p->B85!=6 && p->B92!=51)
	{
        irregular_parameters(p);
        
        if(p->B92==31)
        {
        amplitudes_irregular(p);

        phases_irregular(p);

        pgc->bcast_double(ei,p->wN,0);
        }
        
        if(p->B92==41)
        {
        amplitudes_focused(p);
        phases_focused(p);
        }
	}
    
    if(p->B92==51)
    {
    recon_read(p,pgc);
    recon_parameters(p,pgc);
    }
    
	if(p->B85==4 || p->B85==5 || p->B85==6)
	wavepackets_parameters(p);

    // Print components based on spectrum type
    if(p->B85==11)
        print_components_2d(p);
    else
        print_components(p);
    
    
    if(p->mpirank==0)
    {
    cout<<"Wave_Lib: 1st-order irregular waves . ";
    if(p->B92==51)
    cout<<"wave_recon ";
    
    cout<<endl;
    
    cout<<"Hs: "<<p->wHs<<" Tp: "<<p->wTp<<" wp: "<<p->wwp<<" cp: "<<p->wC<<" water depth: "<<wdt<<endl;
    if(p->B92>40 && p->B92<50)
    cout<<"Focused Wave   xF: "<< p->B81_1 << " yF: " << p->B81_3 <<" tF: "<<p->B81_2<<endl;
    }
    
    singamma = sin((p->B105_1)*(PI/180.0));
    cosgamma = cos((p->B105_1)*(PI/180.0));

    p->Darray(sinhkd,p->wN);

    for(n=0;n<p->wN;++n)
    sinhkd[n] = sinh(ki[n]*wdt);
}

wave_lib_irregular_1st::~wave_lib_irregular_1st()
{
}

// U -------------------------------------------------------------
double wave_lib_irregular_1st::wave_u(lexer *p, double x, double y, double z)
{
    vel=0.0;
	
	for(n=0;n<p->wN;++n)
	Ti[n] = ki[n]*(cosbeta[n]*x + sinbeta[n]*y) - wi[n]*(p->wavetime) - ei[n];
	
	for(n=0;n<p->wN;++n)
    vel += wi[n]*Ai[n]* (cosh(ki[n]*(wdt+z))/sinhkd[n] ) *cos(Ti[n]) * cosbeta[n];
    
    if(p->B130==0)
    vel*=cosgamma;
    
    return vel;
}

double wave_lib_irregular_1st::wave_u_space_sin(lexer *p, double x, double y, double z, int n)
{
	T = sin(ki[n]*(cosbeta[n]*x + sinbeta[n]*y));
	
    vel = wi[n]*Ai[n]* (cosh(ki[n]*(wdt+z))/sinh(ki[n]*wdt) ) * T * cosbeta[n];

    if(p->B130==0)
    vel*=cosgamma;
    
    return vel;
}

double wave_lib_irregular_1st::wave_u_space_cos(lexer *p, double x, double y, double z, int qn)
{
	T = cos(ki[qn]*(cosbeta[qn]*x + sinbeta[qn]*y));
	
    vel = wi[qn]*Ai[qn]* (cosh(ki[qn]*(wdt+z))/sinh(ki[qn]*wdt) ) * T * cosbeta[qn];

    if(p->B130==0)
    vel*=cosgamma;
    
    return vel;
}

double wave_lib_irregular_1st::wave_u_time_sin(lexer *p, int n)
{
	T = sin( -wi[n]*p->wavetime - ei[n]);
	
    return T;
}

double wave_lib_irregular_1st::wave_u_time_cos(lexer *p, int n)
{
    T = cos( -wi[n]*p->wavetime - ei[n]);
	
    return T;
}


// V -------------------------------------------------------------
double wave_lib_irregular_1st::wave_v(lexer *p, double x, double y, double z)
{
    vel=0.0;
	
	for(n=0;n<p->wN;++n)
	Ti[n] = ki[n]*(cosbeta[n]*x + sinbeta[n]*y) - wi[n]*(p->wavetime) - ei[n];
	
	
	for(n=0;n<p->wN;++n)
    vel += wi[n]*Ai[n]* (cosh(ki[n]*(wdt+z))/sinhkd[n] ) * cos(Ti[n]) * sinbeta[n];
	
    if(p->B130==0)
    vel*=singamma;
    
    return vel;
}

double wave_lib_irregular_1st::wave_v_space_sin(lexer *p, double x, double y, double z, int n)
{
	T = sin(ki[n]*(cosbeta[n]*x + sinbeta[n]*y));
	
    vel = wi[n]*Ai[n]* (cosh(ki[n]*(wdt+z))/sinh(ki[n]*wdt) ) * T * sinbeta[n];
	
    if(p->B130==0)
    vel*=singamma;
    
    return vel;
}

double wave_lib_irregular_1st::wave_v_space_cos(lexer *p, double x, double y, double z, int n)
{
	T = cos(ki[n]*(cosbeta[n]*x + sinbeta[n]*y));
	
    vel = wi[n]*Ai[n]* (cosh(ki[n]*(wdt+z))/sinh(ki[n]*wdt) ) * T * sinbeta[n];
	
    if(p->B130==0)
    vel*=singamma;
    
    return vel;
}

double wave_lib_irregular_1st::wave_v_time_sin(lexer *p, int n)
{
	T = sin(- wi[n]*(p->wavetime) - ei[n]);
    
    return T;
}

double wave_lib_irregular_1st::wave_v_time_cos(lexer *p, int n)
{
	T = cos(- wi[n]*(p->wavetime) - ei[n]);
    
    return T;
}


// W -------------------------------------------------------------
double wave_lib_irregular_1st::wave_w(lexer *p, double x, double y, double z)
{
    vel=0.0;
	
	for(n=0;n<p->wN;++n)
	Ti[n] = ki[n]*(cosbeta[n]*x + sinbeta[n]*y) - wi[n]*(p->wavetime) - ei[n];

	for(n=0;n<p->wN;++n)
    vel += wi[n]*Ai[n]* (sinh(ki[n]*(wdt+z))/sinhkd[n]) * sin(Ti[n]);
	
    return vel;
}

double wave_lib_irregular_1st::wave_w_space_sin(lexer *p, double x, double y, double z, int n)
{
	T = sin(ki[n]*(cosbeta[n]*x + sinbeta[n]*y));

    vel = wi[n]*Ai[n]* (sinh(ki[n]*(wdt+z))/sinh(ki[n]*wdt)) * T;
	
    return vel;
}

double wave_lib_irregular_1st::wave_w_space_cos(lexer *p, double x, double y, double z, int n)
{
	T = cos(ki[n]*(cosbeta[n]*x + sinbeta[n]*y));

    vel = wi[n]*Ai[n]* (sinh(ki[n]*(wdt+z))/sinh(ki[n]*wdt)) * T;
	
    return vel;
}

double wave_lib_irregular_1st::wave_w_time_sin(lexer *p, int n)
{
	T = sin(- wi[n]*(p->wavetime) - ei[n]);

    return T;
}

double wave_lib_irregular_1st::wave_w_time_cos(lexer *p, int n)
{
	T = cos(- wi[n]*(p->wavetime) - ei[n]);

    return T;
}

// ETA -------------------------------------------------------------
double wave_lib_irregular_1st::wave_eta(lexer *p, double x, double y)
{
    // Local loop index/accumulator (the static increment::n and member eta
    // had to be spilled around every libm call) and a SIMD-able term loop:
    // with -fno-math-errno -fopenmp-simd GCC maps cos() to libmvec. The
    // components are still summed in the original order.
    const int N = p->wN;
    const double t = p->wavetime;
    double *const __restrict T = Ti;
    
    #pragma omp simd
	for(int q=0;q<N;++q)
	T[q] = Ai[q]*cos(ki[q]*(cosbeta[q]*x + sinbeta[q]*y) - wi[q]*t - ei[q]);
    
    double acc=0.0;
	for(int q=0;q<N;++q)
    acc += T[q];
	
    eta = acc;
    return eta;
}

double wave_lib_irregular_1st::wave_eta_space_sin(lexer *p, double x, double y, int n)
{
	T = sin(ki[n]*(cosbeta[n]*x + sinbeta[n]*y));

    eta =  Ai[n]*T;
	
    return eta;
}

double wave_lib_irregular_1st::wave_eta_space_cos(lexer *p, double x, double y, int n)
{
	T = cos(ki[n]*(cosbeta[n]*x + sinbeta[n]*y));

    eta =  Ai[n]*T;
	
    return eta;
}

double wave_lib_irregular_1st::wave_eta_time_sin(lexer *p, int n)
{
	T = sin(-wi[n]*(p->wavetime) - ei[n]);
	
    return T;
}

double wave_lib_irregular_1st::wave_eta_time_cos(lexer *p, int n)
{
	T = cos(-wi[n]*(p->wavetime) - ei[n]);
	
    return T;
}

// FI -------------------------------------------------------------
double wave_lib_irregular_1st::wave_fi(lexer *p, double x, double y, double z)
{
    const int N = p->wN;
    const double t = p->wavetime;
    const double zz = wdt+z;
    double *const __restrict T = Ti;
    
    #pragma omp simd
    for(int q=0;q<N;++q)
    T[q] = ((wi[q]*Ai[q])/ki[q])*(cosh(ki[q]*zz)/sinhkd[q] ) * sin(ki[q]*(cosbeta[q]*x + sinbeta[q]*y) - wi[q]*t - ei[q]);
    
    double acc=0.0;
    for(int q=0;q<N;++q)
    acc += T[q];
    
    fi = acc;
    return fi;
}
    

double wave_lib_irregular_1st::wave_fi_space_sin(lexer *p, double x, double y, double z, int n)
{
	T = sin(ki[n]*(cosbeta[n]*x + sinbeta[n]*y));
	
    fi = ((wi[n]*Ai[n])/ki[n])*(cosh(ki[n]*(wdt+z))/sinh(ki[n]*wdt) ) * T;
    
    return fi;
}

double wave_lib_irregular_1st::wave_fi_space_cos(lexer *p, double x, double y, double z, int n)
{
    T = cos(ki[n]*(cosbeta[n]*x + sinbeta[n]*y));
	
    fi = ((wi[n]*Ai[n])/ki[n])*(cosh(ki[n]*(wdt+z))/sinh(ki[n]*wdt) ) * T;
    
    return fi;
}

double wave_lib_irregular_1st::wave_fi_time_sin(lexer *p, int n)
{
    T = sin( -wi[n]*p->wavetime - ei[n]);
	
    return T;
}

double wave_lib_irregular_1st::wave_fi_time_cos(lexer *p, int n)
{
    T = cos( -wi[n]*p->wavetime - ei[n]);
	
    return T;
}

void wave_lib_irregular_1st::parameters(lexer *p, ghostcell *pgc)
{

}

void wave_lib_irregular_1st::wave_prestep(lexer *p, ghostcell *pgc)
{
}

// ---- cached-point evaluation ------------------------------------------------
//
// phase_n(x,y,t) = a_n(x,y) + b_n(t),  a_n = k_n (cosbeta_n x + sinbeta_n y),
//                                       b_n = -w_n t - e_n
// cos/sin(a_n) are stored per registered cell, cos/sin(b_n) per step, and
// the sums use the angle-addition formulae. The depth functions are formed
// from one exp per cell and component:
//   cosh(k(d+z))/sinh(kd) = (e^{kz} + e^{-2kd} e^{-kz}) / (1 - e^{-2kd})
//   sinh(k(d+z))/sinh(kd) = (e^{kz} - e^{-2kd} e^{-kz}) / (1 - e^{-2kd})
// which is also better behaved for large kd than cosh()/sinh() directly.
// Results agree with the plain functions to round-off (not bit-identical).

void wave_lib_irregular_1st::wave_cache_points(lexer *p, const std::vector<double> &x, const std::vector<double> &y)
{
    cache_x=x;
    cache_y=y;

    const int N=int(x.size());
    const int M=p->wN;

    cS.assign(size_t(N)*M,0.0);
    sS.assign(size_t(N)*M,0.0);

    for(int q=0;q<N;++q)
    for(int m=0;m<M;++m)
    {
        const double a = ki[m]*(cosbeta[m]*x[q] + sinbeta[m]*y[q]);
        cS[size_t(q)*M+m] = cos(a);
        sS[size_t(q)*M+m] = sin(a);
    }

    cT.assign(M,0.0); sT.assign(M,0.0); ezb.assign(M,0.0);
    em2kd.assign(M,0.0); invden.assign(M,0.0);
    Aeta.assign(M,0.0); Afi.assign(M,0.0); Au.assign(M,0.0); Av.assign(M,0.0); Aw.assign(M,0.0);

    for(int m=0;m<M;++m)
    {
        em2kd[m]  = exp(-2.0*ki[m]*wdt);
        invden[m] = 1.0/(1.0 - em2kd[m]);

        Aeta[m] = Ai[m];
        Afi[m]  = (wi[m]*Ai[m])/ki[m];
        Au[m]   = wi[m]*Ai[m]*cosbeta[m];
        Av[m]   = wi[m]*Ai[m]*sinbeta[m];
        Aw[m]   = wi[m]*Ai[m];
    }

    cache_t=-1.0e300;
}

void wave_lib_irregular_1st::cache_time(lexer *p)
{
    if(p->wavetime==cache_t)
    return;

    const int M=p->wN;
    const double t=p->wavetime;

    for(int m=0;m<M;++m)
    {
        const double b = -wi[m]*t - ei[m];
        cT[m]=cos(b);
        sT[m]=sin(b);
    }

    cache_t=t;
}

double wave_lib_irregular_1st::wave_eta_c(lexer *p, int q)
{
    cache_time(p);

    const int M=p->wN;
    const double *cs=&cS[size_t(q)*M], *ss=&sS[size_t(q)*M];

    double acc=0.0;
    for(int m=0;m<M;++m)
    acc += Aeta[m]*(cs[m]*cT[m] - ss[m]*sT[m]);      // cos(a+b)

    return acc;
}

double wave_lib_irregular_1st::wave_fi_c(lexer *p, int q, double z)
{
    cache_time(p);

    const int M=p->wN;
    const double *cs=&cS[size_t(q)*M], *ss=&sS[size_t(q)*M];
    double *ez=&ezb[0];

    #pragma omp simd
    for(int m=0;m<M;++m)
    ez[m]=exp(ki[m]*z);

    double acc=0.0;
    for(int m=0;m<M;++m)
    {
        const double chr = (ez[m] + em2kd[m]/ez[m])*invden[m];       // cosh ratio
        acc += Afi[m]*chr*(ss[m]*cT[m] + cs[m]*sT[m]);                // sin(a+b)
    }

    return acc;
}

double wave_lib_irregular_1st::wave_u_c(lexer *p, int q, double z)
{
    cache_time(p);

    const int M=p->wN;
    const double *cs=&cS[size_t(q)*M], *ss=&sS[size_t(q)*M];
    double *ez=&ezb[0];

    #pragma omp simd
    for(int m=0;m<M;++m)
    ez[m]=exp(ki[m]*z);

    double acc=0.0;
    for(int m=0;m<M;++m)
    {
        const double chr = (ez[m] + em2kd[m]/ez[m])*invden[m];
        acc += Au[m]*chr*(cs[m]*cT[m] - ss[m]*sT[m]);
    }

    if(p->B130==0)
    acc*=cosgamma;

    return acc;
}

double wave_lib_irregular_1st::wave_v_c(lexer *p, int q, double z)
{
    cache_time(p);

    const int M=p->wN;
    const double *cs=&cS[size_t(q)*M], *ss=&sS[size_t(q)*M];
    double *ez=&ezb[0];

    #pragma omp simd
    for(int m=0;m<M;++m)
    ez[m]=exp(ki[m]*z);

    double acc=0.0;
    for(int m=0;m<M;++m)
    {
        const double chr = (ez[m] + em2kd[m]/ez[m])*invden[m];
        acc += Av[m]*chr*(cs[m]*cT[m] - ss[m]*sT[m]);
    }

    if(p->B130==0)
    acc*=singamma;

    return acc;
}

double wave_lib_irregular_1st::wave_w_c(lexer *p, int q, double z)
{
    cache_time(p);

    const int M=p->wN;
    const double *cs=&cS[size_t(q)*M], *ss=&sS[size_t(q)*M];
    double *ez=&ezb[0];

    #pragma omp simd
    for(int m=0;m<M;++m)
    ez[m]=exp(ki[m]*z);

    double acc=0.0;
    for(int m=0;m<M;++m)
    {
        const double shr = (ez[m] - em2kd[m]/ez[m])*invden[m];       // sinh ratio
        acc += Aw[m]*shr*(ss[m]*cT[m] + cs[m]*sT[m]);
    }

    return acc;
}

void wave_lib_irregular_1st::wave_uvw_c(lexer *p, int q, double z, double &u, double &v, double &w)
{
    cache_time(p);

    const int M=p->wN;
    const double *cs=&cS[size_t(q)*M], *ss=&sS[size_t(q)*M];
    double *ez=&ezb[0];

    #pragma omp simd
    for(int m=0;m<M;++m)
    ez[m]=exp(ki[m]*z);

    double au=0.0, av=0.0, aw=0.0;
    for(int m=0;m<M;++m)
    {
        const double r   = em2kd[m]/ez[m];
        const double chr = (ez[m] + r)*invden[m];
        const double shr = (ez[m] - r)*invden[m];
        const double cph = cs[m]*cT[m] - ss[m]*sT[m];
        const double sph = ss[m]*cT[m] + cs[m]*sT[m];

        au += Au[m]*chr*cph;
        av += Av[m]*chr*cph;
        aw += Aw[m]*shr*sph;
    }

    if(p->B130==0)
    {
        au*=cosgamma;
        av*=singamma;
    }

    u=au; v=av; w=aw;
}
