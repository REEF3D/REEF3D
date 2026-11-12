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

#ifndef FNPF_PROF_H_
#define FNPF_PROF_H_

#include<mpi.h>

class lexer;
class ghostcell;

// FNPF time step profile: wall time of the parts of a step, summed over the steps.
// The parts of fnpf_RK3::start / fnpf_RK4::start are timed at the top level (no nesting); what is not timed
// there (the RK combinations, the copies of K) is rk_rest. report() prints the mean over
// the ranks and the slowest rank per step, every 100 steps and at the end of the run.
namespace fnpf_prof
{
    enum id
    {
        // driver::loop_fnpf
        wavegen, rk_total, print, timestep,
        // fnpf_RK3::start / fnpf_RK4::start
        fsfbc, fsfdisc, sigma, bc3d, gc3d, laplace, fsfwvel, velcalc,
        relax, damping, wetdry, coastline, breaking, inflow, gc2d, body,
        N
    };
    
    extern double t[N];
    extern double tstep;
    extern int steps;
    
    struct scope
    {
        id k;
        double t0;
        explicit scope(id kk) : k(kk), t0(MPI_Wtime()) {}
        ~scope() {t[k] += MPI_Wtime() - t0;}
    };
    
    void step_end(double steptime);
    void report(lexer*, ghostcell*, bool final);
}

#endif
