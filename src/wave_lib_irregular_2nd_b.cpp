/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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
        pgc->bcast_double(ei,p->wN);
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
    cout<<"Hs: "<<p->wHs<<" Tp: "<<p->wTp<<" wp: "<<p->wwp<<endl;
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
	Ti[n] = ki[n]*(cosbeta[n]*x + sinbeta[n]*y) - wi[n]*(p->simtime) - ei[n];
	
	// 1st-order
	for(n=0;n<p->wN;++n)
    vel += wi[n]*Ai[n]*(cosh(ki[n]*(wdt+z))/sinhkd[n] ) * cos(Ti[n]) * cosbeta[n];
    
    // 2nd-order
    for(n=0;n<p->wN-1;++n)
    for(m=n+1;m<p->wN;++m)
    {
    denom1 = Dplus[n][m]*cosh((ki[n]+ki[m])*wdt);
    denom2 = Dminus[n][m]*cosh((ki[n]-ki[m])*wdt);
    denom1 = fabs(denom1)>1.0e-20?denom1:1.0e20;
    denom2 = fabs(denom2)>1.0e-20?denom2:1.0e20;
    
    vel += (ki[n]+ki[m])*Ai[n]*Ai[m]*((Gplus[n][m]*cosh((ki[n]+ki[m])*(z+wdt)))/denom1)*cos(Ti[n]+Ti[m])*(cosbeta[n]*cosbeta[m] + sinbeta[n]*sinbeta[m])
        +  (ki[n]-ki[m])*Ai[n]*Ai[m]*((Gminus[n][m]*cosh((ki[n]-ki[m])*(z+wdt)))/denom2)*cos(Ti[n]-Ti[m])*(cosbeta[n]*cosbeta[m] - sinbeta[n]*sinbeta[m]);
    }
    
    for(n=0;n<p->wN;++n)
    {
     denom3 = Dplus[n][n]*cosh(2.0*ki[n]*wdt); 
     denom3 = fabs(denom3)>1.0e-20?denom3:1.0e20;
     
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
	Ti[n] = ki[n]*(cosbeta[n]*x + sinbeta[n]*y) - wi[n]*(p->simtime) - ei[n];
	
	// 1st-order
	for(n=0;n<p->wN;++n)
    vel += wi[n]*Ai[n]* (cosh(ki[n]*(wdt+z))/sinhkd[n] ) * cos(Ti[n]) * sinbeta[n];
    
    // 2nd-order
    for(n=0;n<p->wN-1;++n)
    for(m=n+1;m<p->wN;++m)
    {
    denom1 = Dplus[n][m]*cosh((ki[n]+ki[m])*wdt);
    denom2 = Dminus[n][m]*cosh((ki[n]-ki[m])*wdt);
    denom1 = fabs(denom1)>1.0e-20?denom1:1.0e20;
    denom2 = fabs(denom2)>1.0e-20?denom2:1.0e20;
    
    vel += (ki[n]+ki[m])*Ai[n]*Ai[m]*((Gplus[n][m]*cosh((ki[n]+ki[m])*(z+wdt)))/denom1)*cos(Ti[n]+Ti[m])*(sinbeta[n]*cosbeta[m] + cosbeta[n]*sinbeta[m])
        +  (ki[n]-ki[m])*Ai[n]*Ai[m]*((Gminus[n][m]*cosh((ki[n]-ki[m])*(z+wdt)))/denom2)*cos(Ti[n]-Ti[m])*(sinbeta[n]*cosbeta[m] - cosbeta[n]*sinbeta[m]);
    }
    
    for(n=0;n<p->wN;++n)
    {
     denom3 = Dplus[n][n]*cosh(2.0*ki[n]*wdt); 
     denom3 = fabs(denom3)>1.0e-20?denom3:1.0e20;
     
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
    Ti[n] = ki[n]*(cosbeta[n]*x + sinbeta[n]*y) - wi[n]*(p->simtime) - ei[n];
    
    // 1st-order
	for(n=0;n<p->wN;++n)
    vel += wi[n]*Ai[n]* (sinh(ki[n]*(wdt+z))/sinhkd[n] ) * sin(Ti[n]);
    
    // 2nd-order
    for(n=0;n<p->wN-1;++n)
    for(m=n+1;m<p->wN;++m)
    {
    denom1 = Dplus[n][m]*cosh((ki[n]+ki[m])*wdt);
    denom2 = Dminus[n][m]*cosh((ki[n]-ki[m])*wdt);
    denom1 = fabs(denom1)>1.0e-20?denom1:1.0e20;
    denom2 = fabs(denom2)>1.0e-20?denom2:1.0e20;
    
    vel += (ki[n]+ki[m])*Ai[n]*Ai[m]*((Gplus[n][m]*sinh((ki[n]+ki[m])*(z+wdt)))/denom1)*sin(Ti[n]+Ti[m])
        +  (ki[n]-ki[m])*Ai[n]*Ai[m]*((Gminus[n][m]*sinh((ki[n]-ki[m])*(z+wdt)))/denom2)*sin(Ti[n]-Ti[m]);
    }
    
    for(n=0;n<p->wN;++n)
    {
     denom3 = Dplus[n][n]*cosh(2.0*ki[n]*wdt); 
     denom3 = fabs(denom3)>1.0e-20?denom3:1.0e20;
     
     vel += ki[n]*Ai[n]*Ai[n]*((Gplus[n][n]*sinh(2.0*ki[n]*(z+wdt)))/denom3)*sin(2.0*Ti[n]);
    }
	
    return vel;
}

double wave_lib_irregular_2nd_b::wave_eta(lexer *p, double x, double y)
{
    eta=0.0;
	
	for(n=0;n<p->wN;++n)
	Ti[n] = ki[n]*(cosbeta[n]*x + sinbeta[n]*y) - wi[n]*(p->simtime) - ei[n];

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
	Ti[n] = ki[n]*(cosbeta[n]*x + sinbeta[n]*y) - wi[n]*(p->simtime) - ei[n];
	
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
     denom3 = pow(sinh(ki[n]*wdt),4.0); 
     denom3 = fabs(denom3)>1.0e-20?denom3:1.0e20;
     
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

