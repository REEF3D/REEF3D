/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

        pgc->bcast_double(ei,p->wN);
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

    print_components(p);
    
    
    if(p->mpirank==0)
    {
    cout<<"Wave Tank: 1st-order irregular waves . ";
    if(p->B92==51)
    cout<<"wave_recon ";
    
    cout<<";  Hs: "<<p->wHs<<" Tp: "<<p->wTp<<" wp: "<<p->wwp<<endl;
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
	Ti[n] = ki[n]*(cosbeta[n]*x + sinbeta[n]*y) - wi[n]*(p->simtime) - ei[n];
	
	for(n=0;n<p->wN;++n)
    if(z<=Ai[n]*cos(Ti[n]))
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
	T = sin( -wi[n]*p->simtime - ei[n]);
	
    return T;
}

double wave_lib_irregular_1st::wave_u_time_cos(lexer *p, int n)
{
    T = cos( -wi[n]*p->simtime - ei[n]);
	
    return T;
}


// V -------------------------------------------------------------
double wave_lib_irregular_1st::wave_v(lexer *p, double x, double y, double z)
{
    vel=0.0;
	
	for(n=0;n<p->wN;++n)
	Ti[n] = ki[n]*(cosbeta[n]*x + sinbeta[n]*y) - wi[n]*(p->simtime) - ei[n];
	
	
	for(n=0;n<p->wN;++n)
    if(z<=Ai[n]*cos(Ti[n]))
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
	T = sin(- wi[n]*(p->simtime) - ei[n]);
    
    return T;
}

double wave_lib_irregular_1st::wave_v_time_cos(lexer *p, int n)
{
	T = cos(- wi[n]*(p->simtime) - ei[n]);
    
    return T;
}


// W -------------------------------------------------------------
double wave_lib_irregular_1st::wave_w(lexer *p, double x, double y, double z)
{
    vel=0.0;
	
	for(n=0;n<p->wN;++n)
	Ti[n] = ki[n]*(cosbeta[n]*x + sinbeta[n]*y) - wi[n]*(p->simtime) - ei[n];

	for(n=0;n<p->wN;++n)
    if(z<=Ai[n]*cos(Ti[n]))
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
	T = sin(- wi[n]*(p->simtime) - ei[n]);

    return T;
}

double wave_lib_irregular_1st::wave_w_time_cos(lexer *p, int n)
{
	T = cos(- wi[n]*(p->simtime) - ei[n]);

    return T;
}

// ETA -------------------------------------------------------------
double wave_lib_irregular_1st::wave_eta(lexer *p, double x, double y)
{
    eta=0.0;

	for(n=0;n<p->wN;++n)
	Ti[n] = ki[n]*(cosbeta[n]*x + sinbeta[n]*y) - wi[n]*(p->simtime) - ei[n];

	for(n=0;n<p->wN;++n)
    eta +=  Ai[n]*cos(Ti[n]);
	
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
	T = sin(-wi[n]*(p->simtime) - ei[n]);
	
    return T;
}

double wave_lib_irregular_1st::wave_eta_time_cos(lexer *p, int n)
{
	T = cos(-wi[n]*(p->simtime) - ei[n]);
	
    return T;
}

// FI -------------------------------------------------------------
double wave_lib_irregular_1st::wave_fi(lexer *p, double x, double y, double z)
{
    fi=0.0;
    
    for(n=0;n<p->wN;++n)
	Ti[n] = ki[n]*(cosbeta[n]*x + sinbeta[n]*y) - wi[n]*(p->simtime) - ei[n];

    for(n=0;n<p->wN;++n)
    fi += ((wi[n]*Ai[n])/ki[n])*(cosh(ki[n]*(wdt+z))/sinhkd[n] ) * sin(Ti[n]);

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
    T = sin( -wi[n]*p->simtime - ei[n]);
	
    return T;
}

double wave_lib_irregular_1st::wave_fi_time_cos(lexer *p, int n)
{
    T = cos( -wi[n]*p->simtime - ei[n]);
	
    return T;
}

void wave_lib_irregular_1st::parameters(lexer *p, ghostcell *pgc)
{

}

void wave_lib_irregular_1st::wave_prestep(lexer *p, ghostcell *pgc)
{
}
