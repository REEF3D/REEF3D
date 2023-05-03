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

#include"wave_lib_Stokes_5th.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

wave_lib_Stokes_5th::wave_lib_Stokes_5th(lexer *p, ghostcell *pgc) : wave_lib_parameters(p,pgc)
{   
    parameters(p,pgc);
    
    if(p->mpirank==0)
    {
    cout<<"Wave Tank: 5th-order Stokes waves; ";
    cout<<"wk: "<<wk<<" ww: "<<ww<<" wf: "<<wf<<" wT: "<<wT<<" wL: "<<wL<<" wdt: "<<wdt<<" kd: "<<wdt*wk<<endl;
    }
    
    singamma = sin((p->B105_1)*(PI/180.0));
    cosgamma = cos((p->B105_1)*(PI/180.0));
}

wave_lib_Stokes_5th::~wave_lib_Stokes_5th()
{
}

// U -------------------------------------------------------------
double wave_lib_Stokes_5th::wave_u(lexer *p, double x, double y, double z)
{
	
	vel = wave_horzvel(p,x,y,z);

    return cosgamma*vel;
}

double wave_lib_Stokes_5th::wave_u_space_sin(lexer *p, double x, double y, double z, int n)
{

	vel = wave_horzvel_space_sin(p,x,y,z,n);

    return cosgamma*vel;
}

double wave_lib_Stokes_5th::wave_u_space_cos(lexer *p, double x, double y, double z, int n)
{

	vel = wave_horzvel_space_cos(p,x,y,z,n);

    return cosgamma*vel;
}

double wave_lib_Stokes_5th::wave_u_time_sin(lexer *p, int n)
{

	vel = wave_horzvel_time_sin(p,n);

    return vel;
}

double wave_lib_Stokes_5th::wave_u_time_cos(lexer *p, int n)
{

	vel = wave_horzvel_time_cos(p,n);

    return vel;
}

// V -------------------------------------------------------------
double wave_lib_Stokes_5th::wave_v(lexer *p, double x, double y, double z)
{

	vel = wave_horzvel(p,x,y,z);

    return singamma*vel;
}

double wave_lib_Stokes_5th::wave_v_space_sin(lexer *p, double x, double y, double z, int n)
{

	vel = wave_horzvel_space_sin(p,x,y,z,n);

    return singamma*vel;
}

double wave_lib_Stokes_5th::wave_v_space_cos(lexer *p, double x, double y, double z, int n)
{

	vel = wave_horzvel_space_cos(p,x,y,z,n);

    return singamma*vel;
}

double wave_lib_Stokes_5th::wave_v_time_sin(lexer *p, int n)
{
	vel = wave_horzvel_time_sin(p,n);

    return vel;
}

double wave_lib_Stokes_5th::wave_v_time_cos(lexer *p, int n)
{
	
	vel = wave_horzvel_time_cos(p,n);

    return vel;
}


// HORZVEL -------------------------------------------------------------
double wave_lib_Stokes_5th::wave_horzvel(lexer *p, double x, double y, double z)
{
	T = wk*x-ww*(p->simtime) + pshift;

    vel = c0*sqrt(9.81/wk)
         *((eps*a11 + pow(eps,3.0)*a31 + pow(eps,5.0)*a51)*cosh(wk*(wdt+z))*cos(T)
         + 2.0*(pow(eps,2.0)*a22 + pow(eps,4.0)*a42)*cosh(2.0*wk*(wdt+z))*cos(2.0*T)
         + 3.0*(pow(eps,3.0)*a33 + pow(eps,5.0)*a53)*cosh(3.0*wk*(wdt+z))*cos(3.0*T)
         + 4.0*(pow(eps,4.0)*a44)*cosh(4.0*wk*(wdt+z))*cos(4.0*T)
         + 5.0*(pow(eps,5.0)*a55)*cosh(5.0*wk*(wdt+z))*cos(5.0*T));

    return vel;
}

double wave_lib_Stokes_5th::wave_horzvel_space_sin(lexer *p, double x, double y, double z, int n)
{
	T = wk*x;

    switch(n)
    {
        case 0: vel = 1.0*c0*sqrt(9.81/wk)*(eps*a11 + pow(eps,3.0)*a31 + pow(eps,5.0)*a51)*cosh(wk*(wdt+z))*sin(1.0*T);
        break;
             
        case 1: vel = 2.0*c0*sqrt(9.81/wk)*(pow(eps,2.0)*a22 + pow(eps,4.0)*a42)*cosh(2.0*wk*(wdt+z))*sin(2.0*T);
        break;
             
        case 2: vel = 3.0*c0*sqrt(9.81/wk)*(pow(eps,3.0)*a33 + pow(eps,5.0)*a53)*cosh(3.0*wk*(wdt+z))*sin(3.0*T);
        break;
             
        case 3: vel = 4.0*c0*sqrt(9.81/wk)*(pow(eps,4.0)*a44)*cosh(4.0*wk*(wdt+z))*sin(4.0*T);
        break;
             
        case 4: vel = 5.0*c0*sqrt(9.81/wk)*(pow(eps,5.0)*a55)*cosh(5.0*wk*(wdt+z))*sin(5.0*T);
        break;
    }

    return vel;
}

double wave_lib_Stokes_5th::wave_horzvel_space_cos(lexer *p, double x, double y, double z, int n)
{

	T = wk*x;

    switch(n)
    {
        case 0: vel = 1.0*c0*sqrt(9.81/wk)*(eps*a11 + pow(eps,3.0)*a31 + pow(eps,5.0)*a51)*cosh(wk*(wdt+z))*cos(1.0*T);
        break;
             
        case 1: vel = 2.0*c0*sqrt(9.81/wk)*(pow(eps,2.0)*a22 + pow(eps,4.0)*a42)*cosh(2.0*wk*(wdt+z))*cos(2.0*T);
        break;
             
        case 2: vel = 3.0*c0*sqrt(9.81/wk)*(pow(eps,3.0)*a33 + pow(eps,5.0)*a53)*cosh(3.0*wk*(wdt+z))*cos(3.0*T);
        break;
             
        case 3: vel = 4.0*c0*sqrt(9.81/wk)*(pow(eps,4.0)*a44)*cosh(4.0*wk*(wdt+z))*cos(4.0*T);
        break;
             
        case 4: vel = 5.0*c0*sqrt(9.81/wk)*(pow(eps,5.0)*a55)*cosh(5.0*wk*(wdt+z))*cos(5.0*T);
        break;
    }

    return vel;
}

double wave_lib_Stokes_5th::wave_horzvel_time_sin(lexer *p, int n)
{
	T = -ww*(p->simtime) + pshift;

    switch(n)
    {
        case 0: vel = sin(1.0*T);
        break;
             
        case 1: vel = sin(2.0*T);
        break;
             
        case 2: vel = sin(3.0*T);
        break;
             
        case 3: vel = sin(4.0*T);
        break;
             
        case 4: vel = sin(5.0*T);
        break;
    }

    return vel;
}

double wave_lib_Stokes_5th::wave_horzvel_time_cos(lexer *p, int n)
{
    T = -ww*(p->simtime) + pshift;

    switch(n)
    {
        case 0: vel = cos(1.0*T);
        break;
             
        case 1: vel = cos(2.0*T);
        break;
             
        case 2: vel = cos(3.0*T);
        break;
             
        case 3: vel = cos(4.0*T);
        break;
             
        case 4: vel = cos(5.0*T);
        break;
    }

    return vel;
}


// W -------------------------------------------------------------
double wave_lib_Stokes_5th::wave_w(lexer *p, double x, double y, double z)
{
	
	T = wk*x-ww*(p->simtime) + pshift;

    vel = c0*sqrt(9.81/wk)
         *((eps*a11 + pow(eps,3.0)*a31 + pow(eps,5.0)*a51)*sinh(wk*(wdt+z))*sin(T)
         + 2.0*(pow(eps,2.0)*a22 + pow(eps,4.0)*a42)*sinh(2.0*wk*(wdt+z))*sin(2.0*T)
         + 3.0*(pow(eps,3.0)*a33 + pow(eps,5.0)*a53)*sinh(3.0*wk*(wdt+z))*sin(3.0*T)
         + 4.0*(pow(eps,4.0)*a44)*sinh(4.0*wk*(wdt+z))*sin(4.0*T)
         + 5.0*(pow(eps,5.0)*a55)*sinh(5.0*wk*(wdt+z))*sin(5.0*T));

    return vel;
}

double wave_lib_Stokes_5th::wave_w_space_sin(lexer *p, double x, double y, double z, int n)
{
	
	T = wk*x;
    
    switch(n)
    {
        case 0: vel = 1.0*c0*sqrt(9.81/wk)*(eps*a11 + pow(eps,3.0)*a31 + pow(eps,5.0)*a51)*sinh(wk*(wdt+z))*sin(1.0*T);
        break;
         
        case 1: vel = 2.0*c0*sqrt(9.81/wk)*(pow(eps,2.0)*a22 + pow(eps,4.0)*a42)*sinh(2.0*wk*(wdt+z))*sin(2.0*T);
        break;
        
        case 2: vel = 3.0*c0*sqrt(9.81/wk)*(pow(eps,3.0)*a33 + pow(eps,5.0)*a53)*sinh(3.0*wk*(wdt+z))*sin(3.0*T);
        break;
        
        case 3: vel = 4.0*c0*sqrt(9.81/wk)*(pow(eps,4.0)*a44)*sinh(4.0*wk*(wdt+z))*sin(4.0*T);
        break;
        
        case 4: vel = 5.0*c0*sqrt(9.81/wk)*(pow(eps,5.0)*a55)*sinh(5.0*wk*(wdt+z))*sin(5.0*T);
        break;
        
        default: vel = 0.0;
        break;
    }

    return vel;
}

double wave_lib_Stokes_5th::wave_w_space_cos(lexer *p, double x, double y, double z, int n)
{
	T = wk*x;
    
    switch(n)
    {
        case 0: vel = 1.0*c0*sqrt(9.81/wk)*(eps*a11 + pow(eps,3.0)*a31 + pow(eps,5.0)*a51)*sinh(wk*(wdt+z))*cos(1.0*T);
        break;
         
        case 1: vel = 2.0*c0*sqrt(9.81/wk)*(pow(eps,2.0)*a22 + pow(eps,4.0)*a42)*sinh(2.0*wk*(wdt+z))*cos(2.0*T);
        break;
        
        case 2: vel = 3.0*c0*sqrt(9.81/wk)*(pow(eps,3.0)*a33 + pow(eps,5.0)*a53)*sinh(3.0*wk*(wdt+z))*cos(3.0*T);
        break;
        
        case 3: vel = 4.0*c0*sqrt(9.81/wk)*(pow(eps,4.0)*a44)*sinh(4.0*wk*(wdt+z))*cos(4.0*T);
        break;
        
        case 4: vel = 5.0*c0*sqrt(9.81/wk)*(pow(eps,5.0)*a55)*sinh(5.0*wk*(wdt+z))*cos(5.0*T);
        break;
        
        default: vel = 0.0;
        break;
    }

    return vel;
}

double wave_lib_Stokes_5th::wave_w_time_sin(lexer *p, int n)
{
	T = -ww*(p->simtime) + pshift;

    switch(n)
    {
        case 0: vel = sin(1.0*T);
        break;
             
        case 1: vel = sin(2.0*T);
        break;
             
        case 2: vel = sin(3.0*T);
        break;
             
        case 3: vel = sin(4.0*T);
        break;
             
        case 4: vel = sin(5.0*T);
        break;
        
        default: vel = 0.0;
        break;
    }

    return vel;
}

double wave_lib_Stokes_5th::wave_w_time_cos(lexer *p, int n)
{
	T = -ww*(p->simtime) + pshift;

    switch(n)
    {
        case 0: vel = cos(1.0*T);
        break;
             
        case 1: vel = cos(2.0*T);
        break;
             
        case 2: vel = cos(3.0*T);
        break;
             
        case 3: vel = cos(4.0*T);
        break;
             
        case 4: vel = cos(5.0*T);
        break;
        
        default: vel = 0.0;
        break;
    }

    return vel;
}

// ETA -------------------------------------------------------------
double wave_lib_Stokes_5th::wave_eta(lexer *p, double x, double y)
{
    double eta;
	
	T = wk*x-ww*(p->simtime) + pshift;

    eta =  (1.0/wk)*((eps + pow(eps,3.0)*b31 - pow(eps,5.0)*(b53 + b55))*cos(T)
                    + (pow(eps,2.0)*b22 + pow(eps,4.0)*b42)*cos(2.0*T)
                    + (-pow(eps,3.0)*b31 + pow(eps,5.0)*b53)*cos(3.0*T)
                    + pow(eps,4.0)*b44*cos(4.0*T)
                    + pow(eps,5.0)*b55*cos(5.0*T));

    return eta;
}

double wave_lib_Stokes_5th::wave_eta_space_sin(lexer *p, double x, double y, int n)
{
	T = wk*x;
    
    switch(n)
    {
        case 0: eta =  (1.0/wk)*(eps + pow(eps,3.0)*b31 - pow(eps,5.0)*(b53 + b55))*sin(1.0*T);
        break;
        
        case 1: eta =  (1.0/wk)*(pow(eps,2.0)*b22 + pow(eps,4.0)*b42)*sin(2.0*T);
        break;
        
        case 2: eta =  (1.0/wk)*(-pow(eps,3.0)*b31 + pow(eps,5.0)*b53)*sin(3.0*T);
        break;
        
        case 3: eta =  (1.0/wk)*pow(eps,4.0)*b44*sin(4.0*T);
        break;
        
        case 4: eta =  (1.0/wk)*pow(eps,5.0)*b55*sin(5.0*T);
        break;
    }
    
    return eta;
}

double wave_lib_Stokes_5th::wave_eta_space_cos(lexer *p, double x, double y, int n)
{
	T = wk*x;
    
    switch(n)
    {
        case 0: eta =  (1.0/wk)*(eps + pow(eps,3.0)*b31 - pow(eps,5.0)*(b53 + b55))*cos(1.0*T);
        break;
        
        case 1: eta =  (1.0/wk)*(pow(eps,2.0)*b22 + pow(eps,4.0)*b42)*cos(2.0*T);
        break;
        
        case 2: eta =  (1.0/wk)*(-pow(eps,3.0)*b31 + pow(eps,5.0)*b53)*cos(3.0*T);
        break;
        
        case 3: eta =  (1.0/wk)*pow(eps,4.0)*b44*cos(4.0*T);
        break;
        
        case 4: eta =  (1.0/wk)*pow(eps,5.0)*b55*cos(5.0*T);
        break;
    }
    
    return eta;
}

double wave_lib_Stokes_5th::wave_eta_time_sin(lexer *p, int n)
{
	T = -ww*(p->simtime) + pshift;

    switch(n)
    {
        case 0: eta = sin(1.0*T);
        break;
             
        case 1: eta = sin(2.0*T);
        break;
             
        case 2: eta = sin(3.0*T);
        break;
             
        case 3: eta = sin(4.0*T);
        break;
             
        case 4: eta = sin(5.0*T);
        break;
    }

    return eta;
}

double wave_lib_Stokes_5th::wave_eta_time_cos(lexer *p, int n)
{
	T = -ww*(p->simtime) + pshift;

    switch(n)
    {
        case 0: eta = cos(1.0*T);
        break;
             
        case 1: eta = cos(2.0*T);
        break;
             
        case 2: eta = cos(3.0*T);
        break;
             
        case 3: eta = cos(4.0*T);
        break;
             
        case 4: eta = cos(5.0*T);
        break;
    }

    return eta;
}

// FI -------------------------------------------------------------
double wave_lib_Stokes_5th::wave_fi(lexer *p, double x, double y, double z)
{
    double fi;
    
    T = wk*x-ww*(p->simtime) + pshift;

    fi = c0*sqrt(9.81/pow(wk,3.0))
         *((eps*a11 + pow(eps,3.0)*a31 + pow(eps,5.0)*a51)*cosh(wk*(wdt+z))*sin(T)
         + 2.0*(pow(eps,2.0)*a22 + pow(eps,4.0)*a42)*cosh(2.0*wk*(wdt+z))*sin(2.0*T)
         + 3.0*(pow(eps,3.0)*a33 + pow(eps,5.0)*a53)*cosh(3.0*wk*(wdt+z))*sin(3.0*T)
         + 4.0*(pow(eps,4.0)*a44)*cosh(4.0*wk*(wdt+z))*sin(4.0*T)
         + 5.0*(pow(eps,5.0)*a55)*cosh(5.0*wk*(wdt+z))*sin(5.0*T));
        
    return fi;
}

void wave_lib_Stokes_5th::wave_fi_precalc_xy_ini(lexer*,int)
{
    
}

void wave_lib_Stokes_5th::wave_fi_precalc_xy(lexer*,double,double,int)
{
    
}

void wave_lib_Stokes_5th::wave_fi_precalc_n(lexer*)
{
    
}

double wave_lib_Stokes_5th::wave_fi_space_sin(lexer *p, double x, double y, double z, int n)
{
    double fi;
    
    T = wk*x;
    
    switch(n)
    {
        case 0: vel = 1.0*c0*sqrt(9.81/pow(wk,3.0))*(eps*a11 + pow(eps,3.0)*a31 + pow(eps,5.0)*a51)*cosh(wk*(wdt+z))*sin(1.0*T);
        break;
         
        case 1: vel = 2.0*c0*sqrt(9.81/pow(wk,3.0))*(pow(eps,2.0)*a22 + pow(eps,4.0)*a42)*cosh(2.0*wk*(wdt+z))*sin(2.0*T);
        break;
        
        case 2: vel = 3.0*c0*sqrt(9.81/pow(wk,3.0))*(pow(eps,3.0)*a33 + pow(eps,5.0)*a53)*cosh(3.0*wk*(wdt+z))*sin(3.0*T);
        break;
        
        case 3: vel = 4.0*c0*sqrt(9.81/pow(wk,3.0))*(pow(eps,4.0)*a44)*cosh(4.0*wk*(wdt+z))*sin(4.0*T);
        break;
        
        case 4: vel = 5.0*c0*sqrt(9.81/pow(wk,3.0))*(pow(eps,5.0)*a55)*cosh(5.0*wk*(wdt+z))*sin(5.0*T);
        break;
        
        default: fi = 0.0;
        break;
    }

    return fi;
}

double wave_lib_Stokes_5th::wave_fi_space_cos(lexer *p, double x, double y, double z, int n)
{
    double fi;
    
    T = wk*x;
    
    switch(n)
    {
        case 0: vel = 1.0*c0*sqrt(9.81/pow(wk,3.0))*(eps*a11 + pow(eps,3.0)*a31 + pow(eps,5.0)*a51)*cosh(wk*(wdt+z))*cos(1.0*T);
        break;
         
        case 1: vel = 2.0*c0*sqrt(9.81/pow(wk,3.0))*(pow(eps,2.0)*a22 + pow(eps,4.0)*a42)*cosh(2.0*wk*(wdt+z))*cos(2.0*T);
        break;
        
        case 2: vel = 3.0*c0*sqrt(9.81/pow(wk,3.0))*(pow(eps,3.0)*a33 + pow(eps,5.0)*a53)*cosh(3.0*wk*(wdt+z))*cos(3.0*T);
        break;
        
        case 3: vel = 4.0*c0*sqrt(9.81/pow(wk,3.0))*(pow(eps,4.0)*a44)*cosh(4.0*wk*(wdt+z))*cos(4.0*T);
        break;
        
        case 4: vel = 5.0*c0*sqrt(9.81/pow(wk,3.0))*(pow(eps,5.0)*a55)*cosh(5.0*wk*(wdt+z))*cos(5.0*T);
        break;
        
        default: fi = 0.0;
        break;
    }

    return fi;
}

double wave_lib_Stokes_5th::wave_fi_time_sin(lexer *p, int n)
{
    T = -ww*(p->simtime) + pshift;

    switch(n)
    {
        case 0: vel = sin(1.0*T);
        break;
             
        case 1: vel = sin(2.0*T);
        break;
             
        case 2: vel = sin(3.0*T);
        break;
             
        case 3: vel = sin(4.0*T);
        break;
             
        case 4: vel = sin(5.0*T);
        break;
        
        default: vel = 0.0;
        break;
    }

    return vel;
}

double wave_lib_Stokes_5th::wave_fi_time_cos(lexer *p, int n)
{
    T = -ww*(p->simtime) + pshift;

    switch(n)
    {
        case 0: vel = cos(1.0*T);
        break;
             
        case 1: vel = cos(2.0*T);
        break;
             
        case 2: vel = cos(3.0*T);
        break;
             
        case 3: vel = cos(4.0*T);
        break;
             
        case 4: vel = cos(5.0*T);
        break;
        
        default: vel = 0.0;
        break;
    }

    return vel;
}

void wave_lib_Stokes_5th::parameters(lexer *p, ghostcell *pgc)
{
    eps = 0.5*wk*wH;


    S = 1.0/cosh(2*wk*wdt);
    C = 1.0 - S;


    a11 = 1.0/sinh(wk*wdt);

    a22 = 3.0*S*S/(2.0*C*C);

    a31 = (-4.0 - 20.0*S + 10.0*S*S -13.0*S*S*S)/(8.0*sinh(wk*wdt)*C*C*C);

    a33 = (-2.0*S*S + 11*S*S*S)/(8.0*sinh(wk*wdt)*C*C*C);

    a42 = (12.0*S - 14.0*S*S - 264.0*S*S*S - 45.0*pow(S,4.0) - 13.0*pow(S,5.0))/(24.0*pow(C,5.0));

    a44 = (10.0*S*S*S - 174.0*pow(S,4.0) + 291.0*pow(S,5.0) + 278.0*pow(S,6.0))/(48.0*(3.0 + 2.0*S)*pow(C,5.0));

    a51 = (-1184.0 + 32.0*S + 13232.0*S*S + 21712.0*S*S*S + 20940.0*pow(S,4.0) + 12554.0*pow(S,5.0)
           -500.0*pow(S,6.0) - 3341.0*pow(S,7.0) - 670.0*pow(S,8.0))/(64.0*sinh(wk*wdt)*(3.0 + 2.0*S)*(4.0 + S)*pow(C,6.0));

    a53 = (4.0*S + 105.0*pow(S,2.0) + 198.0*pow(S,3.0) - 1376.0*pow(S,4.0) - 1302.0*pow(S,5.0) - 117.0*pow(S,6.0) + 58.0*pow(S,7.0))
          /(32.0*sinh(wk*wdt)*(3.0 + 2.0*S)*pow(C,6.0));

    a55 = (-6.0*S*S*S + 272.0*pow(S,4.0) - 1552.0*pow(S,5.0) + 852.0*pow(S,6.0) + 2029.0*pow(S,7.0) + 430.0*pow(S,8.0))
          /(64.0*sinh(wk*wdt)*(3.0 + 2.0*S)*(4.0 + S)*pow(C,6.0));
          


    b22 = ((cosh(wk*wdt)/sinh(wk*wdt))*(1.0 + 2.0*S))/(2.0*C);

    b31 = (-3.0*(1.0 + 3.0*S + 3.0*S*S + 2.0*S*S*S))/(8.0*C*C*C);

    b42 = ((cosh(wk*wdt)/sinh(wk*wdt))*(6.0 - 26.0*S - 182.0*S*S - 204.0*S*S*S - 25.0*pow(S,4.0) + 26*pow(S,5.0)))
          /(6.0*(3.0 + 2.0*S)*pow(C,4.0));

    b44 = ((cosh(wk*wdt)/sinh(wk*wdt))*(24.0 + 92.0*S + 122.0*S*S + 66.0*S*S*S + 67.0*pow(S,4.0) + 34.0*pow(S,5.0)))
          /(24.0*(3.0 + 2.0*S)*pow(C,4.0));

    b53 = (9.0*(132.0 + 17.0*S - 2216.0*S*S - 5897.0*S*S*S - 6292.0*pow(S,4.0) - 2687.0*pow(S,5.0) + 194.0*pow(S,6.0)
                + 467.0*pow(S,7.0) + 82.0*pow(S,8.0)))/(128.0*(3.0 + 2.0*S)*(4.0 + S)*pow(C,6.0));

    b55 = (5.0*(300.0 + 1579.0*S + 3176.0*S*S + 2949.0*S*S*S + 1188.0*pow(S,4.0) + 675.0*pow(S,5.0) + 1326.0*pow(S,6.0)
                + 827.0*pow(S,7.0) + 130.0*pow(S,8.0)))/(384.0*(3.0 + 2.0*S)*(4.0 + S)*pow(C,6.0));


    c0 = sqrt(tanh(wk*wdt));

    c2 = (c0*(2.0 + 7.0*S*S)/(4.0*C*C));

    c4 = (c0*(4.0 + 32.0*S -116.0*S*S - 400.0*S*S*S - 71.0*pow(S,4.0) + 146.0*pow(S,5.0)))/(32.0*pow(C,5.0));


    e2 = (tanh(wk*wdt)*(2.0 + 2.0*S + 5.0*S*S))/(4.0*pow(1.0 - S,2.0));

    e4 = (tanh(wk*wdt)*(8.0 + 12.0*S - 152.0*S*S - 308.0*pow(S,3.0) - 42.0*pow(S,4.0) + 77.0*pow(S,5.0)))/(32.0*pow(1.0 - S, 5.0));
    
}

void wave_lib_Stokes_5th::wave_prestep(lexer *p, ghostcell *pgc)
{
}
