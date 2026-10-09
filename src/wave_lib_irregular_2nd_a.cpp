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

#include"wave_lib_irregular_2nd_a.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

wave_lib_irregular_2nd_a::wave_lib_irregular_2nd_a(lexer *p, ghostcell *pgc) : wave_lib_parameters(p,pgc) 
{ 
    if(p->B85!=4 && p->B85!=5 && p->B85!=6 && p->B92!=52)
	{
        irregular_parameters(p);
        parameters(p,pgc);
        
        if(p->B92==32)
        {
        amplitudes_irregular(p);
        phases_irregular(p);
        pgc->bcast_double(ei,p->wN,0);
        }
        
        if(p->B92==42)
        {
        amplitudes_focused(p);
        phases_focused(p);
        }
	}
	
    if(p->B92==52)
    {
    recon_read(p,pgc);
    recon_parameters(p,pgc);
    parameters(p,pgc);
    }
    
	if(p->B85==4 || p->B85==5 || p->B85==6)
	{
	wavepackets_parameters(p);
	parameters(p,pgc);
	}
    
    print_components(p);
    
    if(p->mpirank==0)
    {
    cout<<"Wave_Lib: 2nd-order irregular waves A"<<endl;
    
    cout<<"Hs: "<<p->wHs<<" Tp: "<<p->wTp<<" wp: "<<p->wwp<<" cp: "<<p->wC<<endl;
    if(p->B92>40 && p->B92<50)
    cout<<"Focused Wave   xF: "<< p->B81_1 << " yF: " << p->B81_3 <<" tF: "<<p->B81_2<<endl;
    }
    
    singamma = sin((p->B105_1)*(PI/180.0));
    cosgamma = cos((p->B105_1)*(PI/180.0));
}

wave_lib_irregular_2nd_a::~wave_lib_irregular_2nd_a()
{
}

double wave_lib_irregular_2nd_a::wave_u(lexer *p, double x, double y, double z)
{
    vel=0.0;
	
	for(n=0;n<p->wN;++n)
	Ti[n] = ki[n]*(cosbeta[n]*x + sinbeta[n]*y) - wi[n]*(p->wavetime) - ei[n];
	
	 // 1st-order
	for(n=0;n<p->wN;++n)
    vel += wi[n]*Ai[n]* (cosh(ki[n]*(wdt+z))/sinh(ki[n]*wdt) ) * cos(Ti[n]) * cosbeta[n];
    
    // 2nd-order
    for(n=0;n<p->wN-1;++n)
    for(m=n+1;m<p->wN;++m)
    {
    denom1 = D1val[n][m];
    denom2 = D2val[n][m];
    
    vel += (Eval[n][m]*cosh((ki[n]-ki[m])*(wdt+z))*(cosbeta[n]*cosbeta[m] + sinbeta[n]*sinbeta[m])*(ki[n]-ki[m]))
        /   denom1
        
        -(Fval[n][m]*cosh((ki[n]+ki[m])*(wdt+z))*cos(Ti[n]+Ti[m])*(cosbeta[n]*cosbeta[m] - sinbeta[n]*sinbeta[m])*(ki[n]-ki[m]))
        /   denom2;
    }
    
    if(p->B130==0)
    vel*=cosgamma;
	
    return vel;
}

double wave_lib_irregular_2nd_a::wave_v(lexer *p, double x, double y, double z)
{
    vel=0.0;
	
	for(n=0;n<p->wN;++n)
	Ti[n] = ki[n]*(cosbeta[n]*x + sinbeta[n]*y) - wi[n]*(p->wavetime) - ei[n];
	
	 // 1st-order
	for(n=0;n<p->wN;++n)
    vel += wi[n]*Ai[n]* (cosh(ki[n]*(wdt+z))/sinh(ki[n]*wdt) ) * cos(Ti[n]) * sinbeta[n];
    
    // 2nd-order
    for(n=0;n<p->wN-1;++n)
    for(m=n+1;m<p->wN;++m)
    {
    denom1 = D1val[n][m];
    denom2 = D2val[n][m];
    
    vel += (Eval[n][m]*cosh((ki[n]-ki[m])*(wdt+z))*cos(Ti[n]-Ti[m])*(sinbeta[n]*cosbeta[m] - cosbeta[n]*sinbeta[m])*(ki[n]-ki[m]))
        /   denom1
        
        -(Fval[n][m]*cosh((ki[n]+ki[m])*(wdt+z))*cos(Ti[n]+Ti[m])*(sinbeta[n]*cosbeta[m] + cosbeta[n]*sinbeta[m])*(ki[n]-ki[m]))
        /   denom2;
    }
    
    if(p->B130==0)
    vel*=singamma;
	
    return vel;
}

double wave_lib_irregular_2nd_a::wave_horzvel(lexer *p, double x, double y, double z)
{
    double vel=0.0;
    
    return vel;
}

double wave_lib_irregular_2nd_a::wave_w(lexer *p, double x, double y, double z)
{
    vel=0.0;

	for(n=0;n<p->wN;++n)
	Ti[n] = ki[n]*(cosbeta[n]*x + sinbeta[n]*y) - wi[n]*(p->wavetime) - ei[n];
    
     // 1st-order
	for(n=0;n<p->wN;++n)
    vel += wi[n]*Ai[n]* (sinh(ki[n]*(wdt+z))/sinh(ki[n]*wdt)) * sin(Ti[n]);
    
    // 2nd-order
    for(n=0;n<p->wN-1;++n)
    for(m=n+1;m<p->wN;++m)
    {
    denom1 = D1val[n][m];
    denom2 = D2val[n][m];
    
    vel += (Eval[n][m]*sinh((ki[n]-ki[m])*(wdt+z))*sin(Ti[n]-Ti[m])*(ki[n]-ki[m]))
        /   denom1
        
        -(Fval[n][m]*sinh(ki[n]+ki[m])*(wdt+z)*sin(Ti[n]+Ti[m])*(ki[n]-ki[m]))
        /   denom2;
    }
	
    return vel;
}

double wave_lib_irregular_2nd_a::wave_eta(lexer *p, double x, double y)
{
    eta=0.0;
		
	for(n=0;n<p->wN;++n)
	Ti[n] = ki[n]*(cosbeta[n]*x + sinbeta[n]*y) - wi[n]*(p->wavetime) - ei[n];

    // 1st-order
	for(n=0;n<p->wN;++n)
    eta +=  Ai[n]*cos(Ti[n]);
    
    // 2nd-order
    for(n=0;n<p->wN-1;++n)
    for(m=n+1;m<p->wN;++m)
    eta +=  ((Ai[n]*Ai[m])/(2.0*fabs(p->W22))) 
        * (Cval[n][m]*cos(Ti[n]-Ti[m]) - Dval[n][m]*cos(Ti[n]+Ti[m]));
	
    return eta;
}

double wave_lib_irregular_2nd_a::wave_fi(lexer *p, double x, double y, double z)
{
    double fi;
    
    return fi;
}

void wave_lib_irregular_2nd_a::parameters(lexer *p, ghostcell *pgc)
{
    p->Darray(Cval,p->wN,p->wN);
    p->Darray(Dval,p->wN,p->wN);
    p->Darray(Eval,p->wN,p->wN);
    p->Darray(Fval,p->wN,p->wN);   
    
    for(n=0;n<p->wN-1;++n)
    for(m=n+1;m<p->wN;++m)
    {
    Cval[n][m] = wave_C(wi[n],wi[m],ki[n],ki[m]);
    Dval[n][m] = wave_D(wi[n],wi[m],ki[n],ki[m]);
    Eval[n][m] = wave_E(wi[n],wi[m],ki[n],ki[m],Ai[n],Ai[m]);
    Fval[n][m] = wave_F(wi[n],wi[m],ki[n],ki[m],Ai[n],Ai[m]);
    }  
    
    // denominators of the velocity terms, once (they were recomputed per pair and evaluation)
    p->Darray(D1val,p->wN,p->wN);
    p->Darray(D2val,p->wN,p->wN);
    
    for(n=0;n<p->wN-1;++n)
    for(m=n+1;m<p->wN;++m)
    {
    denom1 = (fabs(p->W22)*(ki[n]-ki[m])*sinh((ki[n]-ki[m])*wdt) - pow(wi[n]-wi[m],2.0)*cosh((ki[n]-ki[m])*wdt));
    denom2 = (fabs(p->W22)*(ki[n]+ki[m])*sinh((ki[n]+ki[m])*wdt) - pow(wi[n]+wi[m],2.0)*cosh((ki[n]+ki[m])*wdt));
    denom1 = fabs(denom1)>1.0e-20?denom1:1.0e20;
    denom2 = fabs(denom2)>1.0e-20?denom2:1.0e20;
    D1val[n][m] = denom1;
    D2val[n][m] = denom2;
    }
}

double wave_lib_irregular_2nd_a::wave_C(double w1, double w2, double k1, double k2)
{
    double C,a1,a2,denom;

    a1 = 1.0/tanh(k1*wdt);
    a2 = 1.0/tanh(k2*wdt);
    
    denom = (pow(w1,2.0)*(pow(a1,2.0)-1.0) - 2.0*w1*w2*(a1*a2-1.0) + pow(w2,2.0)*(pow(a2,2.0)-1.0));
    
    denom = fabs(denom)>1.0e-20?denom:1.0e20;
    
    C = ((2.0*w1*w2*(w1-w2)*(1.0 + a1*a2) + pow(w1,3.0)*(pow(a1,2.0)-1.0) - pow(w2,3.0)*(pow(a2,2.0)-1.0))*(w1-w2)*(a1*a2-1.0))
        /denom
        - (pow(w1,2.0)+pow(w2,2.0) - w1*w2*(a1*a2+1.0));
        
    return C;
}

double wave_lib_irregular_2nd_a::wave_D(double w1, double w2, double k1, double k2)
{
    double D,a1,a2,denom;

    a1 = 1.0/tanh(k1*wdt);
    a2 = 1.0/tanh(k2*wdt);
    
    denom = (pow(w1,2.0)*(pow(a1,2.0)-1.0) - 2.0*w1*w2*(a1*a2+1.0) + pow(w2,2.0)*(pow(a2,2.0)-1.0));
    
    denom = fabs(denom)>1.0e-20?denom:1.0e20;
    
    D = ((2.0*w1*w2*(w1+w2)*(a1*a2-1.0) + pow(w1,3.0)*(pow(a1,2.0)-1.0) + pow(w2,3.0)*(pow(a2,2.0)-1.0))*(w1+w2)*(a1*a2+1.0))
        /denom
        - (pow(w1,2.0)+pow(w2,2.0) + w1*w2*(a1*a2-1.0));
    
    return D;
}

double wave_lib_irregular_2nd_a::wave_E(double w1, double w2, double k1, double k2, double An, double Am)
{
    double E,a1,a2;

    a1 = 1.0/tanh(k1*wdt);
    a2 = 1.0/tanh(k2*wdt);
    
    E = -0.5*An*Am*(2.0*w1*w2*(w1-w2)*(1.0+a1*a2) + pow(w1,3.0)*(pow(a1,2.0)-1.0) - pow(w2,3.0)*(pow(a2,2.0)-1.0));
    
    return E;
}

double wave_lib_irregular_2nd_a::wave_F(double w1, double w2, double k1, double k2, double An, double Am)
{
    double F,a1,a2;

    a1 = 1.0/tanh(k1*wdt);
    a2 = 1.0/tanh(k2*wdt);
    
    F = -0.5*An*Am*(2.0*w1*w2*(w1+w2)*(1.0-a1*a2) - pow(w1,3.0)*(pow(a1,2.0)-1.0) - pow(w2,3.0)*(pow(a2,2.0)-1.0));
        
    return F;
}

void wave_lib_irregular_2nd_a::wave_prestep(lexer *p, ghostcell *pgc)
{
}


// ---------------------------------------------------------------------
// cached-point evaluation (wave_lib_irregular_2nd_cache.h): the same terms as
// wave_eta, wave_u, wave_v, wave_w above, pair functions from per-component ones
// ---------------------------------------------------------------------

void wave_lib_irregular_2nd_a::wave_cache_points(lexer *p, const std::vector<double> &x, const std::vector<double> &y)
{
    cache_x=x;
    cache_y=y;
    
    const int M=p->wN;
    cc.points(x,y,M,ki,cosbeta,sinbeta);
    cache_coeffs(p);
}

// coefficients of the cached evaluation, once
void wave_lib_irregular_2nd_a::cache_coeffs(lexer *p)
{
    if(coeffs_on)
    return;
    
    const int M=p->wN;
    coeffs_on = true;
    
    fU.assign(M,0.0); fV.assign(M,0.0); fW.assign(M,0.0);
    
    for(n=0;n<M;++n)
    {
        const double a = wi[n]*Ai[n]/sinh(ki[n]*wdt);
        fU[n] = a*cosbeta[n];
        fV[n] = a*sinbeta[n];
        fW[n] = a;
    }
    
    qU1.clear(); qU2.clear(); qV1.clear(); qV2.clear(); qW1.clear(); qW2.clear(); qE1.clear(); qE2.clear();
    
    for(n=0;n<M-1;++n)
    for(m=n+1;m<M;++m)
    {
        const double dk = ki[n]-ki[m];
        const double e1 = Eval[n][m]*dk/D1val[n][m];
        const double e2 = Fval[n][m]*dk/D2val[n][m];
        
        // u: the difference term has no phase factor in wave_u
        qU1.push_back(e1*(cosbeta[n]*cosbeta[m] + sinbeta[n]*sinbeta[m]));
        qU2.push_back(e2*(cosbeta[n]*cosbeta[m] - sinbeta[n]*sinbeta[m]));
        qV1.push_back(e1*(sinbeta[n]*cosbeta[m] - cosbeta[n]*sinbeta[m]));
        qV2.push_back(e2*(sinbeta[n]*cosbeta[m] + cosbeta[n]*sinbeta[m]));
        // w: the sum term is sinh(k_n+k_m)*(d+z) in wave_w
        qW1.push_back(e1);
        qW2.push_back(e2*sinh(ki[n]+ki[m]));
        
        const double ae = (Ai[n]*Ai[m])/(2.0*fabs(p->W22));
        qE1.push_back(ae*Cval[n][m]);
        qE2.push_back(ae*Dval[n][m]);
    }
}

// eta at one point for the times tv (iowave::timeseries)
void wave_lib_irregular_2nd_a::wave_eta_series(lexer *p, double x, double y, const std::vector<double> &tv, std::vector<double> &ev)
{
    cache_coeffs(p);
    wave_lib_irregular_2nd_cache tc;
    tc.points(std::vector<double>(1,x),std::vector<double>(1,y),p->wN,ki,cosbeta,sinbeta);
    ev.assign(tv.size(),0.0);
    
    for(size_t i=0; i<tv.size(); ++i)
    {
        tc.time(tv[i],wi,ei);
        tc.phases(0);
        ev[i] = eta_pairs(p,tc);
    }
}

double wave_lib_irregular_2nd_a::wave_eta_c(lexer *p, int q)
{
    cc.time(p->wavetime,wi,ei);
    cc.phases(q);
    
    return eta_pairs(p,cc);
}

double wave_lib_irregular_2nd_a::eta_pairs(lexer *p, const wave_lib_irregular_2nd_cache &c)
{
    const int M=p->wN;
    const double *C=c.C.data(), *S=c.S.data();
    
    double e=0.0;
    
    for(n=0;n<M;++n)
    e += Ai[n]*C[n];
    
    int k=0;
    for(n=0;n<M-1;++n)
    for(m=n+1;m<M;++m,++k)
    {
        const double cc_ = C[n]*C[m], ss_ = S[n]*S[m];
        e += qE1[k]*(cc_ + ss_) - qE2[k]*(cc_ - ss_);
    }
    
    return e;
}

void wave_lib_irregular_2nd_a::wave_uvw_c(lexer *p, int q, double z, double &u, double &v, double &w)
{
    const int M=p->wN;
    const double s = wdt + z;
    
    if(!cc.vertical(s))
    {
        u = wave_u(p,cache_x[q],cache_y[q],z);
        v = wave_v(p,cache_x[q],cache_y[q],z);
        w = wave_w(p,cache_x[q],cache_y[q],z);
        return;
    }
    
    cc.time(p->wavetime,wi,ei);
    cc.phases(q);
    const double *C=cc.C.data(), *S=cc.S.data(), *E=cc.E.data(), *R=cc.R.data();
    
    double au=0.0, av=0.0, aw=0.0;
    
    for(n=0;n<M;++n)
    {
        const double ch = 0.5*(E[n]+R[n]);
        const double sh = 0.5*(E[n]-R[n]);
        au += fU[n]*ch*C[n];
        av += fV[n]*ch*C[n];
        aw += fW[n]*sh*S[n];
    }
    
    int k=0;
    for(n=0;n<M-1;++n)
    for(m=n+1;m<M;++m,++k)
    {
        const double cc_ = C[n]*C[m], ss_ = S[n]*S[m], sc_ = S[n]*C[m], cs_ = C[n]*S[m];
        const double em = E[n]*R[m], rm = R[n]*E[m];
        const double chm = 0.5*(em+rm), shm = 0.5*(em-rm);
        const double chp = 0.5*(E[n]*E[m] + R[n]*R[m]);
        const double cplus = cc_ - ss_;
        
        au += qU1[k]*chm - qU2[k]*chp*cplus;
        av += qV1[k]*chm*(cc_ + ss_) - qV2[k]*chp*cplus;
        aw += qW1[k]*shm*(sc_ - cs_) - qW2[k]*s*(sc_ + cs_);
    }
    
    if(p->B130==0)
    {
        au*=cosgamma;
        av*=singamma;
    }
    
    u=au; v=av; w=aw;
}

double wave_lib_irregular_2nd_a::wave_u_c(lexer *p, int q, double z)
{
    double u,v,w;
    wave_uvw_c(p,q,z,u,v,w);
    return u;
}

double wave_lib_irregular_2nd_a::wave_v_c(lexer *p, int q, double z)
{
    double u,v,w;
    wave_uvw_c(p,q,z,u,v,w);
    return v;
}

double wave_lib_irregular_2nd_a::wave_w_c(lexer *p, int q, double z)
{
    double u,v,w;
    wave_uvw_c(p,q,z,u,v,w);
    return w;
}
