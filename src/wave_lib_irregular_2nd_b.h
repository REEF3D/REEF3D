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

#include"wave_lib_precalc.h"
#include"wave_lib_parameters.h"
#include"wave_lib_spectrum.h"
#include"increment.h"

#ifndef WAVE_LIB_IRREGULAR_2ND_B_H_
#define WAVE_LIB_IRREGULAR_2ND_B_H_

using namespace std;

class wave_lib_irregular_2nd_b : public wave_lib_precalc, public wave_lib_parameters, public wave_lib_spectrum,
                               public increment
{
public:
    wave_lib_irregular_2nd_b(lexer*, ghostcell*);
	virtual ~wave_lib_irregular_2nd_b();

    double wave_horzvel(lexer*,double,double,double);
    
    virtual double wave_u(lexer*,double,double,double);
    virtual double wave_v(lexer*,double,double,double);
    virtual double wave_w(lexer*,double,double,double);
    virtual double wave_eta(lexer*,double,double);
    virtual double wave_fi(lexer*,double,double,double);
    
    virtual void parameters(lexer*,ghostcell*);
    virtual void wave_prestep(lexer*,ghostcell*);
    
private: 
    double wave_A_plus(double,double,double,double);
	double wave_A_minus(double,double,double,double);
    double wave_D_plus(double,double,double,double);
	double wave_D_minus(double,double,double,double);
	double wave_G_plus(double,double,double,double);
	double wave_G_minus(double,double,double,double);
	double wave_H_plus(double,double,double,double);
	double wave_H_minus(double,double,double,double);
	double wave_F_plus(double,double,double,double);
	double wave_F_minus(double,double,double,double);
    
    double **Aplus,**Aminus,**Dplus,**Dminus,**Gplus,**Gminus,**Hplus,**Hminus,**Fplus,**Fminus;
    
    double *cosh_kpk,*cosh_kmk,*cosh_2k,*sinh_4kh;
    int m;
    double singamma,cosgamma;
    double T,vel,eta,fi;
    double denom1,denom2,denom3;
    
    
    double *sinhkd;
    
};

#endif
