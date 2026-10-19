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

#ifndef WAVE_LIB_FLAP_DOUBLE_H_
#define WAVE_LIB_FLAP_DOUBLE_H_

#include"wave_lib_precalc.h"
#include"wave_lib_parameters.h"
#include"increment.h"
#include"wave_lib_wavemaker2nd.h"

#include<vector>

using namespace std;

class wave_lib_flap_double final : public wave_lib_precalc, public wave_lib_parameters, public increment
{
public:
    // checkpoints (P 48, I 45): the position in the time series of the wavemaker
    void ckpt_state(std::vector<double> &v) override { v.push_back(double(timecount)); v.push_back(double(timecount_old)); v.push_back(double(timecount_z)); v.push_back(double(wtime_old)); }
    void ckpt_restore(const std::vector<double> &v, size_t &k) override { if(k+4>v.size()) return; timecount = decltype(timecount)(v[k++]); timecount_old = decltype(timecount_old)(v[k++]); timecount_z = decltype(timecount_z)(v[k++]); wtime_old = decltype(wtime_old)(v[k++]); }
    wave_lib_flap_double(lexer*, ghostcell*);
	virtual ~wave_lib_flap_double();
    
    double wave_horzvel(lexer*,double,double,double);
    
    double wave_u(lexer*,double,double,double) override final;
    double wave_v(lexer*,double,double,double) override final;
    double wave_w(lexer*,double,double,double) override final;
    double wave_eta(lexer*,double,double) override final;
    double wave_fi(lexer*,double,double,double) override final;
    
    
    double wave_paddle_Q(lexer*,double) override final;
    void wavemaker_setup(lexer*,ghostcell*,double) override final;

    void parameters(lexer*,ghostcell*) override final;
    void wave_prestep(lexer*,ghostcell*) override final;
    
private:
	void read(lexer*, ghostcell*);
	
	int ptnum;
	
	double **kinematics;
	
	double ts,te,val;
	int timecount,timecount_old,timecount_z;
    double wtime_old;
    double singamma,cosgamma;
    
    wave_lib_wavemaker2nd wm2;
    int setup_done;
};

#endif
