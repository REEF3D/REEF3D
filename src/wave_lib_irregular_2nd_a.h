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

#ifndef WAVE_LIB_IRREGULAR_2ND_A_H_
#define WAVE_LIB_IRREGULAR_2ND_A_H_

#include"wave_lib_precalc.h"
#include"wave_lib_parameters.h"
#include"wave_lib_spectrum.h"
#include"increment.h"

using namespace std;

class wave_lib_irregular_2nd_a : public wave_lib_precalc, public wave_lib_parameters, public wave_lib_spectrum,
                               public increment
{
public:
    wave_lib_irregular_2nd_a(lexer*, ghostcell*);
	virtual ~wave_lib_irregular_2nd_a();
    
    double wave_horzvel(lexer*,double,double,double);
    
    virtual double wave_u(lexer*,double,double,double);
    virtual double wave_v(lexer*,double,double,double);
    virtual double wave_w(lexer*,double,double,double);
    virtual double wave_eta(lexer*,double,double);
    virtual double wave_fi(lexer*,double,double,double);
    
    
    virtual void parameters(lexer*,ghostcell*);
    virtual void wave_prestep(lexer*,ghostcell*);
    
private: 
    double wave_C(double,double,double,double);
    double wave_D(double,double,double,double);
    double wave_E(double,double,double,double,double,double);
    double wave_F(double,double,double,double,double,double);
    
    double **Cval,**Dval,**Eval,**Fval;
    int m;
    double singamma,cosgamma;
    double T,vel,eta,fi;
    double denom1,denom2,denom3;
};

#endif
