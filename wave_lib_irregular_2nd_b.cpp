/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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
    cout<<"Wave Tank: 2nd-order irregular waves B; ";
    cout<<";  Hs: "<<p->wHs<<" Tp: "<<p->wTp<<" wp: "<<p->wwp<<endl;
    if(p->B92>40 && p->B92<50)
    cout<<"Focused Wave   xF: "<<p->B81_1<<" tF: "<<p->B81_2<<endl;
    }
    
    singamma = sin((p->B105_1)*(PI/180.0));
    cosgamma = cos((p->B105_1)*(PI/180.0));
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
    vel += wi[n]*Ai[n]* (cosh(ki[n]*(wd+z))/sinh(ki[n]*wd) ) * cos(Ti[n]) * cosbeta[n];
    
    // 2nd-order
    for(n=0;n<p->wN-1;++n)
    for(m=n+1;m<p->wN;++m)
    {
    denom1 = Dplus[n][m]*cosh((ki[n]+ki[m])*wd);
    denom2 = Dminus[n][m]*cosh((ki[n]-ki[m])*wd);
    denom1 = fabs(denom1)>1.0e-20?denom1:1.0e20;
    denom2 = fabs(denom2)>1.0e-20?denom2:1.0e20;
    
    vel += (ki[n]+ki[m])*Ai[n]*Ai[m]*((Gplus[n][m]*cosh((ki[n]+ki[m])*(z+wd)))/denom1)*cos(Ti[n]+Ti[m])*(cosbeta[n]*cosbeta[m] + sinbeta[n]*sinbeta[m])
        +  (ki[n]-ki[m])*Ai[n]*Ai[m]*((Gminus[n][m]*cosh((ki[n]-ki[m])*(z+wd)))/denom2)*cos(Ti[n]-Ti[m])*(cosbeta[n]*cosbeta[m] - sinbeta[n]*sinbeta[m]);
    }
    
    for(n=0;n<p->wN;++n)
    {
     denom3 = Dplus[n][n]*cosh(2.0*ki[n]*wd); 
     denom3 = fabs(denom3)>1.0e-20?denom3:1.0e20;
     
     vel += ki[n]*Ai[n]*Ai[n]*((Gplus[n][n]*cosh(2.0*ki[n]*(z+wd)))/denom3)*cos(2.0*Ti[n]);
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
    vel += wi[n]*Ai[n]* (cosh(ki[n]*(wd+z))/sinh(ki[n]*wd) ) * cos(Ti[n]) * sinbeta[n];
    
    // 2nd-order
    for(n=0;n<p->wN-1;++n)
    for(m=n+1;m<p->wN;++m)
    {
    denom1 = Dplus[n][m]*cosh((ki[n]+ki[m])*wd);
    denom2 = Dminus[n][m]*cosh((ki[n]-ki[m])*wd);
    denom1 = fabs(denom1)>1.0e-20?denom1:1.0e20;
    denom2 = fabs(denom2)>1.0e-20?denom2:1.0e20;
    
    vel += (ki[n]+ki[m])*Ai[n]*Ai[m]*((Gplus[n][m]*cosh((ki[n]+ki[m])*(z+wd)))/denom1)*cos(Ti[n]+Ti[m])*(sinbeta[n]*cosbeta[m] + cosbeta[n]*sinbeta[m])
        +  (ki[n]-ki[m])*Ai[n]*Ai[m]*((Gminus[n][m]*cosh((ki[n]-ki[m])*(z+wd)))/denom2)*cos(Ti[n]-Ti[m])*(sinbeta[n]*cosbeta[m] - cosbeta[n]*sinbeta[m]);
    }
    
    for(n=0;n<p->wN;++n)
    {
     denom3 = Dplus[n][n]*cosh(2.0*ki[n]*wd); 
     denom3 = fabs(denom3)>1.0e-20?denom3:1.0e20;
     
     vel += ki[n]*Ai[n]*Ai[n]*((Gplus[n][n]*cosh(2.0*ki[n]*(z+wd)))/denom3)*cos(2.0*Ti[n]);
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
    vel += wi[n]*Ai[n]* (sinh(ki[n]*(wd+z))/sinh(ki[n]*wd) ) * sin(Ti[n]);
    
    // 2nd-order
    for(n=0;n<p->wN-1;++n)
    for(m=n+1;m<p->wN;++m)
    {
    denom1 = Dplus[n][m]*cosh((ki[n]+ki[m])*wd);
    denom2 = Dminus[n][m]*cosh((ki[n]-ki[m])*wd);
    denom1 = fabs(denom1)>1.0e-20?denom1:1.0e20;
    denom2 = fabs(denom2)>1.0e-20?denom2:1.0e20;
    
    vel += (ki[n]+ki[m])*Ai[n]*Ai[m]*((Gplus[n][m]*sinh((ki[n]+ki[m])*(z+wd)))/denom1)*sin(Ti[n]+Ti[m])
        +  (ki[n]-ki[m])*Ai[n]*Ai[m]*((Gminus[n][m]*sinh((ki[n]-ki[m])*(z+wd)))/denom2)*sin(Ti[n]-Ti[m]);
    }
    
    for(n=0;n<p->wN;++n)
    {
     denom3 = Dplus[n][n]*cosh(2.0*ki[n]*wd); 
     denom3 = fabs(denom3)>1.0e-20?denom3:1.0e20;
     
     vel += ki[n]*Ai[n]*Ai[n]*((Gplus[n][n]*sinh(2.0*ki[n]*(z+wd)))/denom3)*sin(2.0*Ti[n]);
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
    double fi;
    
    return fi;
}

void wave_lib_irregular_2nd_b::parameters(lexer *p, ghostcell *pgc)
{
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
    Dplus[n][m] = wave_D_plus(wi[n],wi[m],ki[n],ki[m]);
    Dminus[n][m] = wave_D_minus(wi[n],wi[m],ki[n],ki[m]);
	Gplus[n][m] = wave_G_plus(wi[n],wi[m],ki[n],ki[m]);
    Gminus[n][m] = wave_G_minus(wi[n],wi[m],ki[n],ki[m]);
	Fplus[n][m] = wave_F_plus(wi[n],wi[m],ki[n],ki[m]);
    Fminus[n][m] = wave_F_minus(wi[n],wi[m],ki[n],ki[m]);
    Hplus[n][m] = wave_H_plus(wi[n],wi[m],ki[n],ki[m]);
    Hminus[n][m] = wave_H_minus(wi[n],wi[m],ki[n],ki[m]);
    }
}


double wave_lib_irregular_2nd_b::wave_D_plus(double w1, double w2, double k1, double k2)
{
    double D;
    
    D = 9.81*(k1+k2)*tanh((k1+k2)*wd) - pow(w1+w2,2.0);
        
    return D;
}

double wave_lib_irregular_2nd_b::wave_D_minus(double w1, double w2, double k1, double k2)
{
    double D;
    
    D = 9.81*(k1-k2)*tanh((k1-k2)*wd) - pow(w1-w2,2.0);
        
    return D;
}

double wave_lib_irregular_2nd_b::wave_G_plus(double w1, double w2, double k1, double k2)
{
	double G,denom1,denom2;
	
	denom1 = 2.0*w1*pow(cosh(k1*wd),2.0);
    denom2 = 2.0*w2*pow(cosh(k2*wd),2.0);
    denom1 = fabs(denom1)>1.0e-20?denom1:1.0e20;
    denom2 = fabs(denom2)>1.0e-20?denom2:1.0e20;
	
	G = -pow(9.81,2.0)*(((k1*k2)/(w1*w2))*(w1+w2)*(1.0-tanh(k1*wd)*tanh(k2*wd)) + (pow(k1,2.0)/denom1 + pow(k2,2.0)/denom2));
	
	return G;	
}

double wave_lib_irregular_2nd_b::wave_G_minus(double w1, double w2, double k1, double k2)
{
	double G,denom1,denom2;
	
	denom1 = 2.0*w1*pow(cosh(k1*wd),2.0);
    denom2 = 2.0*w2*pow(cosh(k2*wd),2.0);
    denom1 = fabs(denom1)>1.0e-20?denom1:1.0e20;
    denom2 = fabs(denom2)>1.0e-20?denom2:1.0e20;
	
	G = -pow(9.81,2.0)*(((k1*k2)/(w1*w2))*(w1-w2)*(1.0+tanh(k1*wd)*tanh(k2*wd)) + (pow(k1,2.0)/denom1 - pow(k2,2.0)/denom2));
	
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
    
    denom1 = (cosh(k1*wd)*cosh(k2*wd));
    denom1 = fabs(denom1)>1.0e-20?denom1:1.0e20;
    
    F = -0.5*9.81*((k1*k2)/(w1*w2))*((pow(cosh((k1-k2)*wd),2.0))/denom1)
        
        + 0.5*(k1*tanh(k1*wd) + k2*tanh(k2*wd)); 
	
	return F;	
}

double wave_lib_irregular_2nd_b::wave_F_minus(double w1, double w2, double k1, double k2)
{
	double F,denom1;
    
    denom1 = (cosh(k1*wd)*cosh(k2*wd));
    denom1 = fabs(denom1)>1.0e-20?denom1:1.0e20;
	
    F = -0.5*9.81*((k1*k2)/(w1*w2))*((pow(cosh((k1-k2)*wd),2.0))/denom1)
        
        + 0.5*(k1*tanh(k1*wd) + k2*tanh(k2*wd)); 
        	
	return F;	
}


