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

#ifndef WAVE_LIB_LINEAR_H_
#define WAVE_LIB_LINEAR_H_

#include<cmath>
#include"wave_lib_precalc.h"
#include"wave_lib_parameters.h"
#include"increment.h"

using namespace std;

class wave_lib_linear final : public wave_lib_precalc, public wave_lib_parameters, public increment
{
public:
    int wave_lexer_fields() const override {return 0;}
    wave_lib_linear(lexer*, ghostcell*);
	virtual ~wave_lib_linear();
    
    double wave_horzvel(lexer*,double,double,double);
    
    double wave_u(lexer*,double,double,double) override final;
    double wave_v(lexer*,double,double,double) override final;
    double wave_w(lexer*,double,double,double) override final;
    double wave_eta(lexer*,double,double) override final;
    double wave_fi(lexer*,double,double,double) override final;
    
    // cached-point evaluation: u, v share the horizontal velocity
    void wave_uvw_c(lexer*, int, double, double&, double&, double&) override final;
    
    void parameters(lexer*,ghostcell*) override final;
    void wave_prestep(lexer*,ghostcell*) override final;
    
    int wave_ncomp() const override final {return 1;}
    void wave_comp(int, double&, double&, double&, double&, double&) const override final;
    double wave_depth0() const override final {return wdt;}
    double wave_depth() const override final {return wdk;}
    void wave_comp_set(int, double, double, double) override final;
    void wave_depth_set(double h) override final {wdk = h;}
    
private:
    double singamma,cosgamma;
    double wsig,wdk;    // intrinsic frequency and depth of the orbital velocities (ww and wdt unless iowave B 530 changes them)
    double wa0,waf;     // amplitude as constructed and the B 530 factor (blocked wave: 0)
    
    // sinh(k d) of the orbital velocities, recomputed when k or d change (B 530)
    double shk_k=-1.0, shk_d=-1.0, shk=0.0;
    double sinhkd()
    {
        if(wk!=shk_k || wdk!=shk_d)
        {
        shk_k = wk;
        shk_d = wdk;
        shk = sinh(wk*wdk);
        }
        return shk;
    }
};

#endif
