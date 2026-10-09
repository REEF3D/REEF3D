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

#ifndef WAVE_LIB_IRREGULAR_1ST_H_
#define WAVE_LIB_IRREGULAR_1ST_H_

#include"wave_lib.h"
#include"wave_lib_parameters.h"
#include"wave_lib_spectrum.h"
#include"increment.h"

using namespace std;

class wave_lib_irregular_1st final : public wave_lib, public wave_lib_parameters, public wave_lib_spectrum,
                               public increment
{
public:
    int wave_lexer_fields() const override {return 0;}
    wave_lib_irregular_1st(lexer*, ghostcell*);
	virtual ~wave_lib_irregular_1st();
    
    double wave_horzvel(lexer*,double,double,double);
    
    // cached-point evaluation (see wave_lib.h): spatial phases per cell and
    // component precomputed, time phases per step, cosh/sinh via one exp
    void wave_cache_points(lexer*, const std::vector<double>&, const std::vector<double>&) override final;
    double wave_eta_c(lexer*, int) override final;
    double wave_fi_c(lexer*, int, double) override final;
    double wave_u_c(lexer*, int, double) override final;
    double wave_v_c(lexer*, int, double) override final;
    double wave_w_c(lexer*, int, double) override final;
    void wave_uvw_c(lexer*, int, double, double&, double&, double&) override final;

    double wave_u(lexer*,double,double,double) override final;
    double wave_u_space_sin(lexer*,double,double,double,int) override final;
    double wave_u_space_cos(lexer*,double,double,double,int) override final;
    double wave_u_time_sin(lexer*,int) override final;
    double wave_u_time_cos(lexer*,int) override final;
    
    double wave_v(lexer*,double,double,double) override final;
    double wave_v_space_sin(lexer*,double,double,double,int) override final;
    double wave_v_space_cos(lexer*,double,double,double,int) override final;
    double wave_v_time_sin(lexer*,int) override final;
    double wave_v_time_cos(lexer*,int) override final;
    
    double wave_w(lexer*,double,double,double) override final;
    double wave_w_space_sin(lexer*,double,double,double,int) override final;
    double wave_w_space_cos(lexer*,double,double,double,int) override final;
    double wave_w_time_sin(lexer*,int) override final;
    double wave_w_time_cos(lexer*,int) override final;
    
    double wave_eta(lexer*,double,double) override final;
    double wave_eta_space_sin(lexer*,double,double,int) override final;
    double wave_eta_space_cos(lexer*,double,double,int) override final;
    double wave_eta_time_sin(lexer*,int) override final;
    double wave_eta_time_cos(lexer*,int) override final;
    
    double wave_fi(lexer*,double,double,double) override final;
    double wave_fi_space_sin(lexer*,double,double,double,int) override final;
    double wave_fi_space_cos(lexer*,double,double,double,int) override final;
    double wave_fi_time_sin(lexer*,int) override final;
    double wave_fi_time_cos(lexer*,int) override final;
    
    
    void parameters(lexer*,ghostcell*) override final;
    void wave_prestep(lexer*,ghostcell*) override final;
    
    // waves on a background (iowave B 530, see wave_lib.h): irregular waves from a spectrum (B 92 31)
    int wave_ncomp() const override final {return bg_ok ? bg_n : 0;}
    void wave_comp(int, double&, double&, double&, double&, double&) const override final;
    double wave_depth0() const override final {return wdt;}
    double wave_depth() const override final {return bg_h;}
    void wave_comp_set(int, double, double, double) override final;
    void wave_depth_set(double) override final;
    void wave_comp_update() override final;
    
private:
    double singamma,cosgamma;    
    double T,vel,eta,fi;
    
    double *sinhkd;

    // cached-point data
    void cache_time(lexer*);
    std::vector<double> cS, sS;                 // cos/sin of the spatial phase, [q*wN+n]
    std::vector<double> cT, sT;                 // cos/sin of the time phase, [n]
    std::vector<double> em2kd, invden;          // exp(-2 k d), 1/(1-exp(-2 k d))
    std::vector<double> Aeta, Afi, Au, Av, Aw;  // component amplitudes of each quantity
    // distinct wave numbers: components of the same frequency in a directional sea share
    // the vertical profile, so exp(k z) is evaluated once per distinct k (results bitwise equal)
    void group_k();
    void vertical_c(double z);
    std::vector<int> iu;                        // component -> distinct k, [n]
    std::vector<double> uK, uem, uinv;          // distinct k, its exp(-2 k d) and 1/(1-exp(-2 k d))
    std::vector<double> uch, ush;               // cosh / sinh ratios at the current z, [distinct k]
    double cache_t=-1.0e300;
    std::vector<double> ph_c, ph_s;             // phases of column ph_q at time ph_t (wave_uvw_c)
    int ph_q=-1;
    double ph_t=-1.0e300;
    double **fixy,*fin;
    
    // B 530 state: the orbital velocities use sigma through sinhkd = sinh(k h_eff) omega/sigma
    // and the cached amplitudes; Ai = Ai0 * amplitude factor
    bool bg_ok=false, bg_dirty=false, bg_used=false;
    int bg_n=0;
    double bg_h=0.0;
    std::vector<double> bg_Ai0, bg_si, bg_af;
    
    int Nw=0, B130v=0;   // p->wN, p->B130 at construction
};

#endif
