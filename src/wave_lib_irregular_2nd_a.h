/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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
#include"wave_lib_irregular_2nd_cache.h"

using namespace std;

class wave_lib_irregular_2nd_a final : public wave_lib_precalc, public wave_lib_parameters, public wave_lib_spectrum,
                               public increment
{
public:
    int wave_lexer_fields() const override {return 1;}
    wave_lib_irregular_2nd_a(lexer*, ghostcell*);
	virtual ~wave_lib_irregular_2nd_a();
    
    
    double wave_u(lexer*,double,double,double) override final;
    double wave_v(lexer*,double,double,double) override final;
    double wave_w(lexer*,double,double,double) override final;
    double wave_eta(lexer*,double,double) override final;
    double wave_fi(lexer*,double,double,double) override final;
    
    
    void parameters(lexer*,ghostcell*) override final;
    void wave_prestep(lexer*,ghostcell*) override final;
    
    // cached-point evaluation (wave_lib_irregular_2nd_cache.h); wave_fi_c stays the direct evaluation
    void wave_cache_points(lexer*, const std::vector<double>&, const std::vector<double>&) override final;
    double wave_eta_c(lexer*, int) override final;
    double wave_u_c(lexer*, int, double) override final;
    double wave_v_c(lexer*, int, double) override final;
    double wave_w_c(lexer*, int, double) override final;
    void wave_uvw_c(lexer*, int, double, double&, double&, double&) override final;
    void wave_eta_series(lexer*, double, double, const std::vector<double>&, std::vector<double>&) override final;
    
private: 
    double wave_C(double,double,double,double);
    double wave_D(double,double,double,double);
    double wave_E(double,double,double,double,double,double);
    double wave_F(double,double,double,double,double,double);
    
    double **Cval,**Dval,**Eval,**Fval;
    double **D1val,**D2val;   // velocity denominators per pair
    int m;
    double singamma,cosgamma;
    double T,vel,eta,fi;
    double denom1,denom2,denom3;
    
    wave_lib_irregular_2nd_cache cc;
    bool coeffs_on=false;
    void cache_coeffs(lexer*);
    double eta_pairs(lexer*, const wave_lib_irregular_2nd_cache&);
    std::vector<double> qU1,qU2,qV1,qV2,qW1,qW2,qE1,qE2;   // pair coefficients, n < m in loop order
    std::vector<double> fU,fV,fW;                          // 1st-order coefficients
};

#endif
