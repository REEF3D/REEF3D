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

    virtual void wave_uvw_c(lexer *p, int q, double z, double &u, double &v, double &w)
    {
        u=wave_u_c(p,q,z);
        v=wave_v_c(p,q,z);
        w=wave_w_c(p,q,z);
    }

    virtual ~wave_lib() = default;

protected:
    std::vector<double> cache_x, cache_y;
};

#endif
