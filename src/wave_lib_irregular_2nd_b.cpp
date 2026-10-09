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

#include"wave_lib_irregular_2nd_b.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

wave_lib_irregular_2nd_b::wave_lib_irregular_2nd_b(lexer *p, ghostcell *pgc) : wave_lib_parameters(p,pgc) 
{ 
    if(p->B85!=4 && p->B85!=5 && p->B85!=6 && p->B92!=53)
	{
        irregular_parameters(p);
        parameters(p,pgc);
        
        if(p->B92==33)
        {
        amplitudes_irregular(p);
        phases_irregular(p);
        pgc->bcast_double(ei,p->wN,0);
        }
        
        if(p->B92==43)
        {
        amplitudes_focused(p);
        phases_focused(p);
        }
	}
    
    if(p->B92==53)
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
    cout<<"Wave_Lib: 2nd-order irregular waves B"<<endl;
    cout<<"Hs: "<<p->wHs<<" Tp: "<<p->wTp<<" wp: "<<p->wwp<<" cp: "<<p->wC<<endl;
    if(p->B92>40 && p->B92<50)
    cout<<"Focused Wave   xF: "<<p->B81_1 << " yF: " << p->B81_3 <<" tF: "<<p->B81_2<<endl;
    }
    
    singamma = sin((p->B105_1)*(PI/180.0));
    cosgamma = cos((p->B105_1)*(PI/180.0));
    
    p->Darray(sinhkd,p->wN);

    for(n=0;n<p->wN;++n)
    sinhkd[n] = sinh(ki[n]*wdt);

}

wave_lib_irregular_2nd_b::~wave_lib_irregular_2nd_b()
{
}

// U -------------------------------------------------------------
double wave_lib_irregular_2nd_b::wave_u(lexer *p, double x, double y, double z)
{
    
    vel=0.0;
	
	for(n=0;n<p->wN;++n)
	Ti[n] = ki[n]*(cosbeta[n]*x + sinbeta[n]*y) - wi[n]*(p->wavetime) - ei[n];
	
	// 1st-order
	for(n=0;n<p->wN;++n)
    vel += wi[n]*Ai[n]*(cosh(ki[n]*(wdt+z))/sinhkd[n] ) * cos(Ti[n]) * cosbeta[n];
    
    // 2nd-order
    for(n=0;n<p->wN-1;++n)
    for(m=n+1;m<p->wN;++m)
    {
    denom1 = P1val[n][m];
    denom2 = P2val[n][m];
    
    vel += (ki[n]+ki[m])*Ai[n]*Ai[m]*((Gplus[n][m]*cosh((ki[n]+ki[m])*(z+wdt)))/denom1)*cos(Ti[n]+Ti[m])*(cosbeta[n]*cosbeta[m] + sinbeta[n]*sinbeta[m])
        +  (ki[n]-ki[m])*Ai[n]*Ai[m]*((Gminus[n][m]*cosh((ki[n]-ki[m])*(z+wdt)))/denom2)*cos(Ti[n]-Ti[m])*(cosbeta[n]*cosbeta[m] - sinbeta[n]*sinbeta[m]);
    }
    
    for(n=0;n<p->wN;++n)
    {
     denom3 = P3val[n];
     
     vel += ki[n]*Ai[n]*Ai[n]*((Gplus[n][n]*cosh(2.0*ki[n]*(z+wdt)))/denom3)*cos(2.0*Ti[n]);
    } 
   
    if(p->B130==0)
    vel*=cosgamma;
	
    return vel;
}


// V -------------------------------------------------------------
double wave_lib_irregular_2nd_b::wave_v(lexer *p, double x, double y, double z)
{
    vel=0.0;
	
	for(n=0;n<p->wN;++n)
	Ti[n] = ki[n]*(cosbeta[n]*x + sinbeta[n]*y) - wi[n]*(p->wavetime) - ei[n];
	
	// 1st-order
	for(n=0;n<p->wN;++n)
    vel += wi[n]*Ai[n]* (cosh(ki[n]*(wdt+z))/sinhkd[n] ) * cos(Ti[n]) * sinbeta[n];
    
    // 2nd-order
    for(n=0;n<p->wN-1;++n)
    for(m=n+1;m<p->wN;++m)
    {
    denom1 = P1val[n][m];
    denom2 = P2val[n][m];
    
    vel += (ki[n]+ki[m])*Ai[n]*Ai[m]*((Gplus[n][m]*cosh((ki[n]+ki[m])*(z+wdt)))/denom1)*cos(Ti[n]+Ti[m])*(sinbeta[n]*cosbeta[m] + cosbeta[n]*sinbeta[m])
        +  (ki[n]-ki[m])*Ai[n]*Ai[m]*((Gminus[n][m]*cosh((ki[n]-ki[m])*(z+wdt)))/denom2)*cos(Ti[n]-Ti[m])*(sinbeta[n]*cosbeta[m] - cosbeta[n]*sinbeta[m]);
    }
    
    for(n=0;n<p->wN;++n)
    {
     denom3 = P3val[n];
     
     vel += ki[n]*Ai[n]*Ai[n]*((Gplus[n][n]*cosh(2.0*ki[n]*(z+wdt)))/denom3)*cos(2.0*Ti[n]);
    }
    
    if(p->B130==0)
    vel*=singamma;
	
    return vel;
}

double wave_lib_irregular_2nd_b::wave_horzvel(lexer *p, double x, double y, double z)
{
    vel=0.0;
    
	
    return vel;
}

double wave_lib_irregular_2nd_b::wave_w(lexer *p, double x, double y, double z)
{
    vel=0.0;

	for(n=0;n<p->wN;++n)
    Ti[n] = ki[n]*(cosbeta[n]*x + sinbeta[n]*y) - wi[n]*(p->wavetime) - ei[n];
    
    // 1st-order
	for(n=0;n<p->wN;++n)
    vel += wi[n]*Ai[n]* (sinh(ki[n]*(wdt+z))/sinhkd[n] ) * sin(Ti[n]);
    
    // 2nd-order
    for(n=0;n<p->wN-1;++n)
    for(m=n+1;m<p->wN;++m)
    {
    denom1 = P1val[n][m];
    denom2 = P2val[n][m];
    
    vel += (ki[n]+ki[m])*Ai[n]*Ai[m]*((Gplus[n][m]*sinh((ki[n]+ki[m])*(z+wdt)))/denom1)*sin(Ti[n]+Ti[m])
        +  (ki[n]-ki[m])*Ai[n]*Ai[m]*((Gminus[n][m]*sinh((ki[n]-ki[m])*(z+wdt)))/denom2)*sin(Ti[n]-Ti[m]);
    }
    
    for(n=0;n<p->wN;++n)
    {
     denom3 = P3val[n];
     
     vel += ki[n]*Ai[n]*Ai[n]*((Gplus[n][n]*sinh(2.0*ki[n]*(z+wdt)))/denom3)*sin(2.0*Ti[n]);
    }
	
    return vel;
}

double wave_lib_irregular_2nd_b::wave_eta(lexer *p, double x, double y)
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
    eta +=  Ai[n]*Ai[m]*Hplus[n][m]*cos(Ti[n]+Ti[m])
          + Ai[n]*Ai[m]*Hminus[n][m]*cos(Ti[n]-Ti[m]);
    
    for(n=0;n<p->wN;++n)
    eta +=  Ai[n]*Ai[n]*Hplus[n][n]*cos(2.0*Ti[n]);
	
    return eta;
}

double wave_lib_irregular_2nd_b::wave_fi(lexer *p, double x, double y, double z)
{
    double fi=0.0;
    
    for(n=0;n<p->wN;++n)
	Ti[n] = ki[n]*(cosbeta[n]*x + sinbeta[n]*y) - wi[n]*(p->wavetime) - ei[n];
	
	// 1st-order
	for(n=0;n<p->wN;++n)
    fi +=  ((wi[n]*Ai[n])/ki[n])*(cosh(ki[n]*(wdt+z))/sinhkd[n] ) * sin(Ti[n]);
    
    
    // 2nd-order
    for(n=0;n<p->wN-1;++n)
    for(m=n+1;m<p->wN;++m)
    {
    fi += Ai[n]*Ai[m]*((Aplus[n][m]*cosh((ki[n]+ki[m])*(z+wdt))))*sin(Ti[n]+Ti[m])*(sinbeta[n]*cosbeta[m] + cosbeta[n]*sinbeta[m])
        + Ai[n]*Ai[m]*((Aminus[n][m]*cosh((ki[n]-ki[m])*(z+wdt))))*sin(Ti[n]-Ti[m])*(sinbeta[n]*cosbeta[m] - cosbeta[n]*sinbeta[m]);
    }
    
    for(n=0;n<p->wN;++n)
    {
     denom3 = P4val[n];
     
     fi += (3.0/8.0)*wi[n]*Ai[n]*Ai[n]*((cosh(2.0*ki[n]*(z+wdt)))/denom3)*sin(2.0*Ti[n]);
    }
    
    return fi;
}

void wave_lib_irregular_2nd_b::parameters(lexer *p, ghostcell *pgc)
{
    p->Darray(Aplus,p->wN,p->wN);
    p->Darray(Aminus,p->wN,p->wN);
    p->Darray(Dplus,p->wN,p->wN);
    p->Darray(Dminus,p->wN,p->wN);
	p->Darray(Gplus,p->wN,p->wN);
    p->Darray(Gminus,p->wN,p->wN);
	p->Darray(Hplus,p->wN,p->wN);
    p->Darray(Hminus,p->wN,p->wN);
	p->Darray(Fplus,p->wN,p->wN);
    p->Darray(Fminus,p->wN,p->wN);
    
    for(n=0;n<p->wN;++n)
    for(m=0;m<p->wN;++m)
    {
    Aplus[n][m] = wave_A_plus(wi[n],wi[m],ki[n],ki[m]);
    Aminus[n][m] = wave_A_minus(wi[n],wi[m],ki[n],ki[m]);
    Dplus[n][m] = wave_D_plus(wi[n],wi[m],ki[n],ki[m]);
    Dminus[n][m] = wave_D_minus(wi[n],wi[m],ki[n],ki[m]);
	Gplus[n][m] = wave_G_plus(wi[n],wi[m],ki[n],ki[m]);
    Gminus[n][m] = wave_G_minus(wi[n],wi[m],ki[n],ki[m]);
	Fplus[n][m] = wave_F_plus(wi[n],wi[m],ki[n],ki[m]);
    Fminus[n][m] = wave_F_minus(wi[n],wi[m],ki[n],ki[m]);
    Hplus[n][m] = wave_H_plus(wi[n],wi[m],ki[n],ki[m]);
    Hminus[n][m] = wave_H_minus(wi[n],wi[m],ki[n],ki[m]);
    
    //cout<<"k: "<<ki[n]<<" "<<ki[m]<<" w: "<<wi[n]<<" "<<wi[m]<<" H+-: "<<Hplus[n][m]<<" "<<Hminus[n][m]<<" F+-: "<<Fplus[n][m]<<" "<<Fminus[n][m]<<endl;
    }
    
    
    // denominators of the 2nd-order terms, once (they were recomputed per pair and evaluation)
    p->Darray(P1val,p->wN,p->wN);
    p->Darray(P2val,p->wN,p->wN);
    p->Darray(P3val,p->wN);
    p->Darray(P4val,p->wN);
    
    for(n=0;n<p->wN-1;++n)
    for(m=n+1;m<p->wN;++m)
    {
    denom1 = Dplus[n][m]*cosh((ki[n]+ki[m])*wdt);
    denom2 = Dminus[n][m]*cosh((ki[n]-ki[m])*wdt);
    denom1 = fabs(denom1)>1.0e-20?denom1:1.0e20;
    denom2 = fabs(denom2)>1.0e-20?denom2:1.0e20;
    P1val[n][m] = denom1;
    P2val[n][m] = denom2;
    }
    
    for(n=0;n<p->wN;++n)
    {
     denom3 = Dplus[n][n]*cosh(2.0*ki[n]*wdt); 
     denom3 = fabs(denom3)>1.0e-20?denom3:1.0e20;
     P3val[n] = denom3;
     denom3 = pow(sinh(ki[n]*wdt),4.0); 
     denom3 = fabs(denom3)>1.0e-20?denom3:1.0e20;
     P4val[n] = denom3;
    }
    
    p->Darray(cosh_kpk,p->wN*p->wN);
    p->Darray(cosh_kmk,p->wN*p->wN);
    p->Darray(cosh_2k,p->wN*p->wN);
    p->Darray(sinh_4kh,p->wN*p->wN);
    
    int count=0;
    for(n=0;n<p->wN-1;++n)
    for(m=n+1;m<p->wN;++m)
    {
        
        +count;
    }
}

double wave_lib_irregular_2nd_b::wave_A_plus(double w1, double w2, double k1, double k2)
{
    double A;
    
    double denom1,denom2;
	
	denom1 = -wave_D_plus(w1,w2,k1,k2);
    denom2 = tanh(k1*wdt)*tanh(k2*wdt);
    
    denom1 = fabs(denom1)>1.0e-20?denom1:1.0e20;
    denom2 = fabs(denom1)>1.0e-20?denom1:1.0e20;
	
	A = -((w1*w2)*(w1+w2)/denom1)*(1.0-1.0/denom2) + (0.5/denom1)*(pow(w1,3.0)/pow(sinh(k1*wdt),2.0) + pow(w2,3.0)/pow(sinh(k2*wdt),2.0));
	
    return A;
}

double wave_lib_irregular_2nd_b::wave_A_minus(double w1, double w2, double k1, double k2)
{
    double A;
    
    double denom1,denom2;
	
	denom1 = -wave_D_minus(w1,w2,k1,k2);
    denom2 = tanh(k1*wdt)*tanh(k2*wdt);
    
    denom1 = fabs(denom1)>1.0e-20?denom1:1.0e20;
    denom2 = fabs(denom1)>1.0e-20?denom1:1.0e20;
	
	A = ((w1*w2)*(w1-w2)/denom1)*(1.0+1.0/denom2) + (0.5/denom1)*(pow(w1,3.0)/pow(sinh(k1*wdt),2.0) - pow(w2,3.0)/pow(sinh(k2*wdt),2.0));
	
    return A;
}

double wave_lib_irregular_2nd_b::wave_D_plus(double w1, double w2, double k1, double k2)
{
    double D;
    
    D = 9.81*(k1+k2)*tanh((k1+k2)*wdt) - pow(w1+w2,2.0);
        
    return D;
}

double wave_lib_irregular_2nd_b::wave_D_minus(double w1, double w2, double k1, double k2)
{
    double D;
    
    D = 9.81*(k1-k2)*tanh((k1-k2)*wdt) - pow(w1-w2,2.0);
        
    return D;
}

double wave_lib_irregular_2nd_b::wave_G_plus(double w1, double w2, double k1, double k2)
{
	double G,denom1,denom2;
	
	denom1 = 2.0*w1*pow(cosh(k1*wdt),2.0);
    denom2 = 2.0*w2*pow(cosh(k2*wdt),2.0);
    denom1 = fabs(denom1)>1.0e-20?denom1:1.0e20;
    denom2 = fabs(denom2)>1.0e-20?denom2:1.0e20;
	
	G = -pow(9.81,2.0)*(((k1*k2)/(w1*w2))*(w1+w2)*(1.0-tanh(k1*wdt)*tanh(k2*wdt)) + (pow(k1,2.0)/denom1 + pow(k2,2.0)/denom2));
	
	return G;	
}

double wave_lib_irregular_2nd_b::wave_G_minus(double w1, double w2, double k1, double k2)
{
	double G,denom1,denom2;
	
	denom1 = 2.0*w1*pow(cosh(k1*wdt),2.0);
    denom2 = 2.0*w2*pow(cosh(k2*wdt),2.0);
    denom1 = fabs(denom1)>1.0e-20?denom1:1.0e20;
    denom2 = fabs(denom2)>1.0e-20?denom2:1.0e20;
	
	G = -pow(9.81,2.0)*(((k1*k2)/(w1*w2))*(w1-w2)*(1.0+tanh(k1*wdt)*tanh(k2*wdt)) + (pow(k1,2.0)/denom1 - pow(k2,2.0)/denom2));
	
	return G;	
}

double wave_lib_irregular_2nd_b::wave_H_plus(double w1, double w2, double k1, double k2)
{
	double H,denom1;
    
    denom1 = wave_D_plus(w1,w2,k1,k2);
    denom1 = fabs(denom1)>1.0e-20?denom1:1.0e20;
	
    H = (w1+w2)*(1.0/9.81)*(wave_G_plus(w1,w2,k1,k2)/denom1) + wave_F_plus(w1,w2,k1,k2);
	
	return H;	
}

double wave_lib_irregular_2nd_b::wave_H_minus(double w1, double w2, double k1, double k2)
{
	double H,denom1;
    
    denom1 = wave_D_minus(w1,w2,k1,k2);
    denom1 = fabs(denom1)>1.0e-20?denom1:1.0e20;
	
    H = (w1-w2)*(1.0/9.81)*(wave_G_minus(w1,w2,k1,k2)/denom1) + wave_F_minus(w1,w2,k1,k2);
	
	return H;	
}

double wave_lib_irregular_2nd_b::wave_F_plus(double w1, double w2, double k1, double k2)
{
	double F,denom1;
    
    denom1 = (cosh(k1*wdt)*cosh(k2*wdt));
    denom1 = fabs(denom1)>1.0e-20?denom1:1.0e20;
    
    F = -0.5*9.81*((k1*k2)/(w1*w2))*((pow(cosh((k1-k2)*wdt),2.0))/denom1)
        
        + 0.5*(k1*tanh(k1*wdt) + k2*tanh(k2*wdt)); 
	
	return F;	
}

double wave_lib_irregular_2nd_b::wave_F_minus(double w1, double w2, double k1, double k2)
{
	double F,denom1;
    
    denom1 = (cosh(k1*wdt)*cosh(k2*wdt));
    denom1 = fabs(denom1)>1.0e-20?denom1:1.0e20;
	
    F = -0.5*9.81*((k1*k2)/(w1*w2))*((pow(cosh((k1-k2)*wdt),2.0))/denom1)
        
        + 0.5*(k1*tanh(k1*wdt) + k2*tanh(k2*wdt)); 
        	
	return F;	
}

void wave_lib_irregular_2nd_b::wave_prestep(lexer *p, ghostcell *pgc)
{
}


// ---------------------------------------------------------------------
// cached-point evaluation (wave_lib_irregular_2nd_cache.h): the same terms as
// wave_eta, wave_fi, wave_u, wave_v, wave_w above, pair functions from per-component ones
// ---------------------------------------------------------------------

void wave_lib_irregular_2nd_b::wave_cache_points(lexer *p, const std::vector<double> &x, const std::vector<double> &y)
{
    cache_x=x;
    cache_y=y;
    
    const int M=p->wN;
    cc.points(x,y,M,ki,cosbeta,sinbeta);
    cache_coeffs(p);
}

// coefficients of the cached evaluation, once
void wave_lib_irregular_2nd_b::cache_coeffs(lexer *p)
{
    if(coeffs_on)
    return;
    
    const int M=p->wN;
    coeffs_on = true;
    
    fU.assign(M,0.0); fV.assign(M,0.0); fW.assign(M,0.0); fF.assign(M,0.0);
    dU.assign(M,0.0); dF.assign(M,0.0); dH.assign(M,0.0);
    
    for(n=0;n<M;++n)
    {
        const double a = wi[n]*Ai[n]/sinhkd[n];
        fU[n] = a*cosbeta[n];
        fV[n] = a*sinbeta[n];
        fW[n] = a;
        fF[n] = ((wi[n]*Ai[n])/ki[n])/sinhkd[n];
        
        dU[n] = ki[n]*Ai[n]*Ai[n]*Gplus[n][n]/P3val[n];
        dF[n] = (3.0/8.0)*wi[n]*Ai[n]*Ai[n]/P4val[n];
        dH[n] = Ai[n]*Ai[n]*Hplus[n][n];
    }
    
    qUp.clear(); qUm.clear(); qVp.clear(); qVm.clear(); qWp.clear(); qWm.clear();
    qFp.clear(); qFm.clear(); qHp.clear(); qHm.clear();
    
    for(n=0;n<M-1;++n)
    for(m=n+1;m<M;++m)
    {
        const double aa = Ai[n]*Ai[m];
        const double gp = (ki[n]+ki[m])*aa*Gplus[n][m]/P1val[n][m];
        const double gm = (ki[n]-ki[m])*aa*Gminus[n][m]/P2val[n][m];
        
        qUp.push_back(gp*(cosbeta[n]*cosbeta[m] + sinbeta[n]*sinbeta[m]));
        qUm.push_back(gm*(cosbeta[n]*cosbeta[m] - sinbeta[n]*sinbeta[m]));
        qVp.push_back(gp*(sinbeta[n]*cosbeta[m] + cosbeta[n]*sinbeta[m]));
        qVm.push_back(gm*(sinbeta[n]*cosbeta[m] - cosbeta[n]*sinbeta[m]));
        qWp.push_back(gp);
        qWm.push_back(gm);
        qFp.push_back(aa*Aplus[n][m]*(sinbeta[n]*cosbeta[m] + cosbeta[n]*sinbeta[m]));
        qFm.push_back(aa*Aminus[n][m]*(sinbeta[n]*cosbeta[m] - cosbeta[n]*sinbeta[m]));
        qHp.push_back(aa*Hplus[n][m]);
        qHm.push_back(aa*Hminus[n][m]);
    }
}

// eta at one point for the times tv (iowave::timeseries)
void wave_lib_irregular_2nd_b::wave_eta_series(lexer *p, double x, double y, const std::vector<double> &tv, std::vector<double> &ev)
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

double wave_lib_irregular_2nd_b::wave_eta_c(lexer *p, int q)
{
    cc.time(p->wavetime,wi,ei);
    cc.phases(q);
    
    return eta_pairs(p,cc);
}

double wave_lib_irregular_2nd_b::eta_pairs(lexer *p, const wave_lib_irregular_2nd_cache &c)
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
        e += qHp[k]*(cc_ - ss_) + qHm[k]*(cc_ + ss_);
    }
    
    for(n=0;n<M;++n)
    e += dH[n]*(C[n]*C[n] - S[n]*S[n]);
    
    return e;
}

double wave_lib_irregular_2nd_b::wave_fi_c(lexer *p, int q, double z)
{
    const int M=p->wN;
    const double s = wdt + z;
    
    if(!cc.vertical(s))
    return wave_fi(p,cache_x[q],cache_y[q],z);
    
    cc.time(p->wavetime,wi,ei);
    cc.phases(q);
    const double *C=cc.C.data(), *S=cc.S.data(), *E=cc.E.data(), *R=cc.R.data();
    
    double f=0.0;
    
    for(n=0;n<M;++n)
    f += fF[n]*0.5*(E[n]+R[n])*S[n];
    
    int k=0;
    for(n=0;n<M-1;++n)
    for(m=n+1;m<M;++m,++k)
    {
        const double sc_ = S[n]*C[m], cs_ = C[n]*S[m];
        const double chm = 0.5*(E[n]*R[m] + R[n]*E[m]);
        const double chp = 0.5*(E[n]*E[m] + R[n]*R[m]);
        f += qFp[k]*chp*(sc_ + cs_) + qFm[k]*chm*(sc_ - cs_);
    }
    
    for(n=0;n<M;++n)
    f += dF[n]*0.5*(E[n]*E[n] + R[n]*R[n])*(2.0*S[n]*C[n]);
    
    return f;
}

void wave_lib_irregular_2nd_b::wave_uvw_c(lexer *p, int q, double z, double &u, double &v, double &w)
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
        const double ep = E[n]*E[m], rp = R[n]*R[m];
        const double chm = 0.5*(em+rm), shm = 0.5*(em-rm);
        const double chp = 0.5*(ep+rp), shp = 0.5*(ep-rp);
        const double cplus = cc_ - ss_, cminus = cc_ + ss_;
        
        au += qUp[k]*chp*cplus + qUm[k]*chm*cminus;
        av += qVp[k]*chp*cplus + qVm[k]*chm*cminus;
        aw += qWp[k]*shp*(sc_ + cs_) + qWm[k]*shm*(sc_ - cs_);
    }
    
    for(n=0;n<M;++n)
    {
        const double e2 = E[n]*E[n], r2 = R[n]*R[n];
        const double c2 = C[n]*C[n] - S[n]*S[n];
        au += dU[n]*0.5*(e2+r2)*c2;
        av += dU[n]*0.5*(e2+r2)*c2;
        aw += dU[n]*0.5*(e2-r2)*(2.0*S[n]*C[n]);
    }
    
    if(p->B130==0)
    {
        au*=cosgamma;
        av*=singamma;
    }
    
    u=au; v=av; w=aw;
}

double wave_lib_irregular_2nd_b::wave_u_c(lexer *p, int q, double z)
{
    double u,v,w;
    wave_uvw_c(p,q,z,u,v,w);
    return u;
}

double wave_lib_irregular_2nd_b::wave_v_c(lexer *p, int q, double z)
{
    double u,v,w;
    wave_uvw_c(p,q,z,u,v,w);
    return v;
}

double wave_lib_irregular_2nd_b::wave_w_c(lexer *p, int q, double z)
{
    double u,v,w;
    wave_uvw_c(p,q,z,u,v,w);
    return w;
}
