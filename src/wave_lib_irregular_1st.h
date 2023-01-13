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

#include"wave_lib.h"
#include"wave_lib_parameters.h"
#include"wave_lib_spectrum.h"
#include"increment.h"

#ifndef WAVE_LIB_IRREGULAR_1ST_H_
#define WAVE_LIB_IRREGULAR_1ST_H_

using namespace std;

class wave_lib_irregular_1st : public wave_lib, public wave_lib_parameters, public wave_lib_spectrum,
                               public increment
{
public:
    wave_lib_irregular_1st(lexer*, ghostcell*);
	virtual ~wave_lib_irregular_1st();
    
    double wave_horzvel(lexer*,double,double,double);
    
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
    virtual double wave_fi_space_sin(lexer*,double,double,double,int);
    virtual double wave_fi_space_cos(lexer*,double,double,double,int);
    virtual double wave_fi_time_sin(lexer*,int);
    virtual double wave_fi_time_cos(lexer*,int);
    
    
    virtual void parameters(lexer*,ghostcell*);
    virtual void wave_prestep(lexer*,ghostcell*);
    
private:
    double singamma,cosgamma;    
    double T,vel,eta,fi;
    
    double *sinhkd;
    double **fixy,*fin;
};

#endif
