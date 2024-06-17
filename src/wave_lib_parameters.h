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

#ifndef WAVE_LIB_PARAMETERS_H_
#define WAVE_LIB_PARAMETERS_H_

class lexer;
class fdm;
class ghostcell;

#include"increment.h"

using namespace std;

class wave_lib_parameters : public increment
{
public:
    wave_lib_parameters(lexer*, ghostcell*);
	virtual ~wave_lib_parameters();
    
    double sinhfunc(double);
    double coshfunc(double);
    
    double sinfunc(double);
    double cosfunc(double);
    
    double teta;
    double wk,ww,wdt,wa,wH,wL,wf,wT,wL0,k0,S0;
    double wk_temp,ww_temp,wL_temp,wT_temp,wf_temp;
    
    
    double eps,c0,c2,c4; 
    double S,C;
    double wC,ubar;
    double wS;
    
    double X0;
	
    const double pshift;
	
    
private: 

	
    int wtype;
    double diff;
    
    double f,r,s;
    int factorial,q;
    const int order;
    
    double *factcos;

};

#endif
