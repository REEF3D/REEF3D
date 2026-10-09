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

#ifndef WAVE_LIB_IRREGULAR_2ND_B_H_
#define WAVE_LIB_IRREGULAR_2ND_B_H_

#include"wave_lib_precalc.h"
#include"wave_lib_parameters.h"
#include"wave_lib_spectrum.h"
#include"increment.h"
#include"wave_lib_irregular_2nd_cache.h"

using namespace std;

class wave_lib_irregular_2nd_b final : public wave_lib_precalc, public wave_lib_parameters, public wave_lib_spectrum,
                               public increment
{
public:
    int wave_lexer_fields() const override {return 0;}
    wave_lib_irregular_2nd_b(lexer*, ghostcell*);
	virtual ~wave_lib_irregular_2nd_b();
    
    double wave_u(lexer*,double,double,double) override final;
    double wave_v(lexer*,double,double,double) override final;
    double wave_w(lexer*,double,double,double) override final;
    double wave_eta(lexer*,double,double) override final;
    double wave_fi(lexer*,double,double,double) override final;
    
    void parameters(lexer*,ghostcell*) override final;
    void wave_prestep(lexer*,ghostcell*) override final;
    
    // cached-point evaluation (wave_lib_irregular_2nd_cache.h)
    void wave_cache_points(lexer*, const std::vector<double>&, const std::vector<double>&) override final;
    double wave_eta_c(lexer*, int) override final;
    double wave_fi_c(lexer*, int, double) override final;
    double wave_u_c(lexer*, int, double) override final;
    double wave_v_c(lexer*, int, double) override final;
    double wave_w_c(lexer*, int, double) override final;
    void wave_uvw_c(lexer*, int, double, double&, double&, double&) override final;
    void wave_eta_series(lexer*, double, double, const std::vector<double>&, std::vector<double>&) override final;
    
private:
    void terms(lexer*);
    void direct(lexer*, double, double);
    void rotate(lexer*, double&, double&);
    
    double singamma,cosgamma;
    
    wave_lib_irregular_2nd_terms tt;      // second-order theory
    wave_lib_irregular_2nd_cache cc, dc;  // cached points / direct evaluation
    
    int Nw=0, B130v=0;   // p->wN, p->B130 at construction
};

#endif
