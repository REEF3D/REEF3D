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

#ifndef WAVE_LIB_STOKES_5TH_H_
#define WAVE_LIB_STOKES_5TH_H_

#include"wave_lib.h"
#include"increment.h"

using namespace std;

class wave_lib_Stokes_5th : public wave_lib, public increment
{
public:
    wave_lib_Stokes_5th(lexer*, ghostcell*);
	virtual ~wave_lib_Stokes_5th();
    
    double wave_horzvel(lexer*,double,double,double);
    virtual double wave_horzvel_space_sin(lexer*,double,double,double,int);
    virtual double wave_horzvel_space_cos(lexer*,double,double,double,int);
    virtual double wave_horzvel_time_sin(lexer*,int);
    virtual double wave_horzvel_time_cos(lexer*,int);
    
    virtual double wave_u(lexer*,double,double,double);
    virtual double wave_u_space_sin(lexer*,double,double,double,int);
    virtual double wave_u_space_cos(lexer*,double,double,double,int);
    virtual double wave_u_time_sin(lexer*,int);
    virtual double wave_u_time_cos(lexer*,int);
    
    virtual double wave_v(lexer*,double,double,double);
    virtual double wave_v_space_sin(lexer*,double,double,double,int);
    virtual double wave_v_space_cos(lexer*,double,double,double,int);
    virtual double wave_v_time_sin(lexer*,int);
    virtual double wave_v_time_cos(lexer*,int);
    
    virtual double wave_w(lexer*,double,double,double);
    virtual double wave_w_space_sin(lexer*,double,double,double,int);
    virtual double wave_w_space_cos(lexer*,double,double,double,int);
    virtual double wave_w_time_sin(lexer*,int);
    virtual double wave_w_time_cos(lexer*,int);
    
    virtual double wave_eta(lexer*,double,double);
    virtual double wave_eta_space_sin(lexer*,double,double,int);
    virtual double wave_eta_space_cos(lexer*,double,double,int);
    virtual double wave_eta_time_sin(lexer*,int);
    virtual double wave_eta_time_cos(lexer*,int);
    
    virtual double wave_fi(lexer*,double,double,double);
    virtual void wave_fi_precalc_xy_ini(lexer*,int);
    virtual void wave_fi_precalc_xy(lexer*,double,double,int);
    virtual void wave_fi_precalc_n(lexer*);
    virtual double wave_fi_space_sin(lexer*,double,double,double,int);
    virtual double wave_fi_space_cos(lexer*,double,double,double,int);
    virtual double wave_fi_time_sin(lexer*,int);
    virtual double wave_fi_time_cos(lexer*,int);
    
    
    virtual void wave_parameters(lexer*,ghostcell*);
    virtual void parameters(lexer*,ghostcell*);
    virtual void wave_prestep(lexer*,ghostcell*);
    
private:
    double a11,a22,a31,a33,a42,a44,a51,a53,a55;
    double b22,b31,b42,b44,b53,b55;
    double e2,e4;
    double singamma,cosgamma;
    double vel,eta,fi,T;
    
    int wtype;
    double diff;
    double teta;
    double wk,ww,wdt,wa,wH,wL,wf,wT,wL0,k0,S0;
    double wk_temp,ww_temp,wL_temp,wT_temp,wT_test,wf_temp;
    
    
    double eps,c0,c2,c4; 
    double S,C;
    double wC,ubar;
    double wS;
    
    double X0;
	
    const double pshift;
};

#endif
