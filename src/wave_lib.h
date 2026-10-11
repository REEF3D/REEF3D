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

#ifndef WAVE_LIB_H_
#define WAVE_LIB_H_

#include<vector>

class lexer;
class fdm;
class ghostcell;
class field;

using namespace std;

class wave_lib
{
public:
    virtual ~wave_lib() = default;

    virtual double wave_u(lexer*,double,double,double)=0;
    virtual double wave_u_space_sin(lexer*,double,double,double,int)=0;
    virtual double wave_u_space_cos(lexer*,double,double,double,int)=0;
    virtual double wave_u_time_sin(lexer*,int)=0;
    virtual double wave_u_time_cos(lexer*,int)=0;
    
    virtual double wave_v(lexer*,double,double,double)=0;
    virtual double wave_v_space_sin(lexer*,double,double,double,int)=0;
    virtual double wave_v_space_cos(lexer*,double,double,double,int)=0;
    virtual double wave_v_time_sin(lexer*,int)=0;
    virtual double wave_v_time_cos(lexer*,int)=0;
    
    virtual double wave_w(lexer*,double,double,double)=0;
    virtual double wave_w_space_sin(lexer*,double,double,double,int)=0;
    virtual double wave_w_space_cos(lexer*,double,double,double,int)=0;
    virtual double wave_w_time_sin(lexer*,int)=0;
    virtual double wave_w_time_cos(lexer*,int)=0;
    
    virtual double wave_eta(lexer*,double,double)=0;
    virtual double wave_eta_space_sin(lexer*,double,double,int)=0;
    virtual double wave_eta_space_cos(lexer*,double,double,int)=0;
    virtual double wave_eta_time_sin(lexer*,int)=0;
    virtual double wave_eta_time_cos(lexer*,int)=0;
    
    virtual double wave_fi(lexer*,double,double,double)=0;
    virtual double wave_fi_space_sin(lexer*,double,double,double,int)=0;
    virtual double wave_fi_space_cos(lexer*,double,double,double,int)=0;
    virtual double wave_fi_time_sin(lexer*,int)=0;
    virtual double wave_fi_time_cos(lexer*,int)=0;
    
    
    // second-order moving-paddle BC terms Q = X_z*phi_z - X*phi_xx at x=0
    // for piston/flap wavemakers (z relative to still water level)
    virtual double wave_paddle_Q(lexer*,double) {return 0.0;}

    virtual void parameters(lexer*,ghostcell*)=0;
    virtual void wave_prestep(lexer*,ghostcell*)=0;

    // ---- cached-point evaluation ------------------------------------------
    // The relaxation-zone cells are fixed, so iowave registers their
    // horizontal coordinates once (wave_cache_points) and afterwards evaluates
    // by cell index. The defaults simply call the plain functions at the
    // stored coordinates, so every wave type works unchanged;
    // wave_lib_irregular_1st overrides them with precomputed spatial phases.
    virtual void wave_cache_points(lexer*, const std::vector<double> &x, const std::vector<double> &y)
    {
        cache_x=x;
        cache_y=y;
    }

    virtual double wave_eta_c(lexer *p, int q) {return wave_eta(p,cache_x[q],cache_y[q]);}
    virtual double wave_fi_c(lexer *p, int q, double z) {return wave_fi(p,cache_x[q],cache_y[q],z);}
    virtual double wave_u_c(lexer *p, int q, double z) {return wave_u(p,cache_x[q],cache_y[q],z);}
    virtual double wave_v_c(lexer *p, int q, double z) {return wave_v(p,cache_x[q],cache_y[q],z);}
    virtual double wave_w_c(lexer *p, int q, double z) {return wave_w(p,cache_x[q],cache_y[q],z);}

    virtual void wave_uvw_c(lexer *p, int q, double z, double &u, double &v, double &w);

    // eta at one point for the times tv (iowave::timeseries, REEF3D_Log-Wave); the default
    // evaluates wave_eta at each time, the 2nd-order theories use their cached evaluation
    virtual void wave_eta_series(lexer *p, double x, double y, const std::vector<double> &tv, std::vector<double> &ev);

    // ---- waves on a background (iowave, B 530) ----------------------------
    // A wave on a tide or current keeps the absolute frequency omega of each
    // component; iowave re-solves k on the background depth h_eff (and with
    // Doppler from omega = sigma + k U_n) and hands the state back. The
    // libraries that support it evaluate eta and the phase with k, the orbital
    // velocities with sigma and sinh(k h_eff); the height above the bed stays
    // wdt + z. A component blocked by an opposing current gets amplitude
    // factor 0.
    //   wave_ncomp        number of components (0: not supported)
    //   wave_comp         component n: k, omega, sigma, direction beta [rad] relative
    //                     to the library's B 105 frame, amplitude factor
    //   wave_depth0/depth construction depth wdt / current h_eff
    //   wave_comp_set, wave_depth_set, then wave_comp_update once to refresh derived data
    virtual int wave_ncomp() const {return 0;}
    virtual void wave_comp(int n, double &k, double &omega, double &sigma, double &beta, double &af) const {}
    virtual double wave_depth0() const {return 0.0;}
    virtual double wave_depth() const {return 0.0;}
    virtual void wave_comp_set(int n, double k, double sigma, double af) {}
    virtual void wave_depth_set(double h) {}
    virtual void wave_comp_update() {}

    // lexer fields the evaluation reads (iowave redesign, step 2): 0 none (the wave parameters are
    // members, set at construction), 1 only wN and B130 (irregular theories), 2 the whole wave
    // context (default). wave_field swaps a source's context into the lexer only as far as needed.
    virtual int wave_lexer_fields() const {return 2;}

protected:
    std::vector<double> cache_x, cache_y;
};

#endif
