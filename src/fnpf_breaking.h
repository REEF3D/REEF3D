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

#ifndef FNPF_BREAKING_H_
#define FNPF_BREAKING_H_

#include"fnpf_fsf.h"
#include"sliceint4.h"
#include"slice4.h"

class fnpf_laplace;
class field;
class fnpf_convection;
class fnpf_ddx;
class fnpf_etadisc;
class fnpf_coastline;
class solver2D;

using namespace std;

class fnpf_breaking : public fnpf_fsf, public increment 
{
public:
	fnpf_breaking(lexer*, fdm_fnpf*, ghostcell*);
	virtual ~fnpf_breaking();
    
    void breaking(lexer*,fdm_fnpf*,ghostcell*,slice&,slice&,slice&,double);
    
    void breaking_baquet(lexer*,fdm_fnpf*,ghostcell*,slice&,slice&,slice&,double);
    
    void breaking_kennedy(lexer*,fdm_fnpf*,ghostcell*,slice&,slice&,slice&,double);
    
    // romero
    void breaking_romero(lexer*,fdm_fnpf*,ghostcell*,slice&,slice&,slice&,double);

    
    
    void filter(lexer*, fdm_fnpf*,ghostcell*, slice&);

    double ivel,jvel,kvel;
    
private:
    // Step 1: Compute local wave steepness S(i,j)
    void compute_steepness(lexer*, fdm_fnpf*, ghostcell*);

    // Step 2: Compute local phase speed c(i,j)
    void compute_phase_speed(lexer*, fdm_fnpf*, ghostcell*);

    // Step 3: Detect breaking onset/cessation
    void detect_breaking(lexer*, fdm_fnpf*, ghostcell*);

    // Step 4: Compute eddy viscosity via energy matching
    void compute_eddy_viscosity(lexer*, fdm_fnpf*, ghostcell*);

    // --- Helper: local wavenumber estimation via zero-crossing ---
    void compute_local_wavenumber(lexer*, fdm_fnpf*, ghostcell*);

    // --- Helper: Hilbert-transform-based envelope and phase ---
    void compute_hilbert_steepness(lexer*, fdm_fnpf*, ghostcell*);
    
    
    double rb3(lexer*,double);
    double rb4(lexer*,double);
    
    double dist3,dist4,db;
    
    double visc;
    
    sliceint4 bx,by;
    
    slice4 eta_t, t_break;
    slice4 B_coeff;
    
    slice4 S_local, b_strength, c_phase, k_local, L_break, B_ramp;
    slice4 eta_old,brk_flag,k_smooth;
    
    int count_n;
    
    double eta_I;       // onset threshold coefficient
    double eta_F;       // cessation threshold coefficient
    double delta_b;     // mixing length coefficient
    double T_star;      // ramp-up duration coefficient
    double alpha_eta;   // KFSBC dissipation scaling (0 = off)
    int ini_done;       // First call flag
    
    
    double a_coeff;       // Romero a coefficient
    double S0;            // onset threshold steepness
    double alpha_L;       // breaking zone width fraction
    double vb_max;      // maximum eddy viscosity limiter
    double eta_t_thresh;  // optional eta_t pre-filter threshold coefficient
    
    int steep_method;     // steepness computation method
    int cspeed_method;    // phase speed computation method
    
};

#endif
