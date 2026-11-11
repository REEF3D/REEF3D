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
Architect: Hans Bihs
--------------------------------------------------------------------*/

#ifndef FNPF_TIMESTEP_H_
#define FNPF_TIMESTEP_H_

#include"timestep.h"
#include"increment.h"
#include<vector>

class fdm_fnpf;
class lexer;
class ghostcell;

using namespace std;

class fnpf_timestep : public increment
{
public:
	fnpf_timestep(lexer*);
	virtual ~fnpf_timestep();
	void start(fdm_fnpf*, lexer*,ghostcell*);
	void ini(fdm_fnpf*, lexer*,ghostcell*);

    // N 48: the rule of the time step, 1 sqrt(g h_max) (RK3, RK4), 2 stability step of the grid, 3 accuracy step of
    // theta-Robin (N 48 1 with A 310 5)
    static int rule(lexer*);
    
    // N 48 1: the gravity wave speed of the CFL criterion (also used by fnpf_amr::timestep)
    static double wave_speed(lexer*, double u, double h, double dx, double hmax);
    
    // N 48 2: frequency of the grid cut-off mode with the discrete Dirichlet-to-Neumann operator of the sigma grid,
    // kh2 the symbol of the horizontal cds2 Laplacian at the Nyquist mode, h the water depth
    static double omega_discrete(lexer*, double kh2, double h);
    static double rk_limit(lexer*);                        // linear stability bound of the RK scheme on the imaginary axis
    static void dtn_table_ini(lexer*);                     // the table g(kH) of the discrete DtN (once, printed on rank 0)
    
    // N 52, N 53: accuracy cap dt <= N52/N53
    static void accuracy_cap(lexer*);
    static double theta_tref(lexer*);                       // theta-Robin, N 48 1: period of the reference wave
    static double theta_error(lexer*, double, double, double); // theta-Robin, N 48 1: phase-speed error of theta-Robin
    static double theta_step(lexer*, double);               // theta-Robin, N 48 1: the step of the accuracy target N 55


private:
    static double dtn_discrete(lexer*, double q);          // dimensionless discrete DtN g(q) = H D, q = kH
    static double dtn_lookup(double q);
    static int sz_weights(lexer*, double*);                // weights of the free-surface sz stencil (A311)
    
    static std::vector<double> tab_lq, tab_g;              // log10(q), g(q)
    static double tab_lqmin, tab_lqmax, tab_dlq, gmax_ratio;
    static int dtn_ok;
    
	double sqd,wallu,wallv,wallw;
	double cu,cv,cw;
	const double epsi;
	double isor,jsor,ksor;
	double irsm,jrsm,krsm;
    const double maxtimestep, c0_orig;
    double dx;

};

#endif
