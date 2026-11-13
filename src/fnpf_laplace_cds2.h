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

#ifndef LAPLACE_FNPF_CDS2_H_
#define LAPLACE_FNPF_CDS2_H_

#include"fnpf_laplace.h"
#include"increment.h"
#include<vector>

class fnpf_bed_update;

using namespace std;

class fnpf_laplace_cds2 final : public fnpf_laplace, public increment
{
public:
    fnpf_laplace_cds2 (lexer*);
	virtual ~fnpf_laplace_cds2();

    void start(lexer *,fdm_fnpf*,ghostcell*,solver*,fnpf_fsf*,double*,slice&) override final;
    
    // mesh refinement (fnpf_amr): assemble M and rhs, but do not solve
    bool assemble_only = false;
    
    // A 430: set during the stage solves of RK3 and RK4 (the tolerance relative to the potential without its
    // mean); the psi solves of the body loads (fnpf_psi_tol) and theta-Robin keep N 44
    static bool stage_tol;
    
    // A 321: the start value of a stage solve of RK3 / RK4, extrapolated in time. kind 1: the first solve of the
    // step, for the time t^n + c dt, from the solutions at t^n (f) and t^n-1: f + c dt/dt_old (f - f^n-1);
    // kind 2: f is the solution at t^n + dt/2, for t^n + dt: 2 f - f^n; kind 3: f is the solution at t^n + dt,
    // for t^n + dt/2: (f + f^n)/2. It is set at the unknowns after the assembly, which keeps the f of the
    // last solve in the explicit sigma cross terms: the solution is the same within the tolerance of the solve.
    // For the next solve only.
    void start_value(lexer*, int kind, double c);
    
    // theta-Robin (A 310 5): if set, the coupling mt of the top rows to the free-surface node (moved to rhs as
    // -mt Fifsf) is stored at their row index; the other rows are not written
    double *mtsurf = nullptr;

    // mesh refinement, A 324 > 0: the bed ghost nodes and the halo of f before the next Picard
    // assembly, as start() does between its own solves
    void picard_bc(lexer*, fdm_fnpf*, ghostcell*, double*, fnpf_fsf*);
    
private:
    fnpf_bed_update *pbed;
    
    double ab,denom;
    double starttime,endtime;
    int gcval;
    
    std::vector<double> fold;
    
    // A 321
    std::vector<double> fn;      // the solution at t^n (t^n-1 until the first solve of the step)
    double sv_a = 1.0, sv_b = 0.0, sv_dt = 0.0;
    bool sv_on = false, sv_store = false;
};

#endif
