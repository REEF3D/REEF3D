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
        pgc->bcast_double(ei,p->wN);
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
	Ti[n] = ki[n]*(cosbeta[n]*x + sinbeta[n]*y) - wi[n]*(p->simtime) - ei[n];
	
	 // 1st-order
	for(n=0;n<p->wN;++n)
    vel += wi[n]*Ai[n]* (cosh(ki[n]*(wdt+z))/sinh(ki[n]*wdt) ) * cos(Ti[n]) * cosbeta[n];
    
    // 2nd-order
    for(n=0;n<p->wN-1;++n)
    for(m=n+1;m<p->wN;++m)
    {
    denom1 = (fabs(p->W22)*(ki[n]-ki[m])*sinh((ki[n]-ki[m])*wdt) - pow(wi[n]-wi[m],2.0)*cosh((ki[n]-ki[m])*wdt));
    denom2 = (fabs(p->W22)*(ki[n]+ki[m])*sinh((ki[n]+ki[m])*wdt) - pow(wi[n]+wi[m],2.0)*cosh((ki[n]+ki[m])*wdt));
    denom1 = fabs(denom1)>1.0e-20?denom1:1.0e20;
    denom2 = fabs(denom2)>1.0e-20?denom2:1.0e20;
    
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
	Ti[n] = ki[n]*(cosbeta[n]*x + sinbeta[n]*y) - wi[n]*(p->simtime) - ei[n];
	
	 // 1st-order
	for(n=0;n<p->wN;++n)
    vel += wi[n]*Ai[n]* (cosh(ki[n]*(wdt+z))/sinh(ki[n]*wdt) ) * cos(Ti[n]) * sinbeta[n];
    
    // 2nd-order
    for(n=0;n<p->wN-1;++n)
    for(m=n+1;m<p->wN;++m)
    {
    denom1 = (fabs(p->W22)*(ki[n]-ki[m])*sinh((ki[n]-ki[m])*wdt) - pow(wi[n]-wi[m],2.0)*cosh((ki[n]-ki[m])*wdt));
    denom2 = (fabs(p->W22)*(ki[n]+ki[m])*sinh((ki[n]+ki[m])*wdt) - pow(wi[n]+wi[m],2.0)*cosh((ki[n]+ki[m])*wdt));
    denom1 = fabs(denom1)>1.0e-20?denom1:1.0e20;
    denom2 = fabs(denom2)>1.0e-20?denom2:1.0e20;
    
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
	Ti[n] = ki[n]*(cosbeta[n]*x + sinbeta[n]*y) - wi[n]*(p->simtime) - ei[n];
    
     // 1st-order
	for(n=0;n<p->wN;++n)
    vel += wi[n]*Ai[n]* (sinh(ki[n]*(wdt+z))/sinh(ki[n]*wdt)) * sin(Ti[n]);
    
    // 2nd-order
    for(n=0;n<p->wN-1;++n)
    for(m=n+1;m<p->wN;++m)
    {
    denom1 = (fabs(p->W22)*(ki[n]-ki[m])*sinh((ki[n]-ki[m])*wdt) - pow(wi[n]-wi[m],2.0)*cosh((ki[n]-ki[m])*wdt));
    denom2 = (fabs(p->W22)*(ki[n]+ki[m])*sinh((ki[n]+ki[m])*wdt) - pow(wi[n]+wi[m],2.0)*cosh((ki[n]+ki[m])*wdt));
    denom1 = fabs(denom1)>1.0e-20?denom1:1.0e20;
    denom2 = fabs(denom2)>1.0e-20?denom2:1.0e20;
    
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
	Ti[n] = ki[n]*(cosbeta[n]*x + sinbeta[n]*y) - wi[n]*(p->simtime) - ei[n];

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

