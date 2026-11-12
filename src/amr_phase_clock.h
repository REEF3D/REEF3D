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

#ifndef AMR_PHASE_CLOCK_H_
#define AMR_PHASE_CLOCK_H_

#include"amr_load.h"
#include<mpi.h>

//  Exclusive clocks of the parts of a solve (the composite projection of CFD AMR): every part is a
//  scope; a scope inside another one takes its time out of the outer one, so the parts add up to
//  the time spent in the outermost scopes.  Per part the wall-clock time and the time in blocking
//  MPI calls within it (amr_clock::mpi, mpi_clock.cpp).  Two MPI_Wtime calls per scope.  The inner
//  scopes (root false) count only inside an outer one: the same routine called outside the solve
//  (a fill of the velocities) is not counted.

template<int N>
struct amr_phase_clock
{
    double t[N] = {};       // wall-clock time per part [s], cumulative
    double m[N] = {};       // of which in blocking MPI calls

    void enter(int k)
    {
        charge();
        if(depth<MAXD)
        stack[depth] = k;
        ++depth;
    }

    bool active() const { return depth>0; }

    void leave()
    {
        charge();
        if(depth>0)
        --depth;
    }

private:
    static const int MAXD = 64;
    int stack[MAXD] = {};
    int depth = 0;
    double tlast = 0.0, mlast = 0.0;

    void charge()
    {
        const double now = MPI_Wtime(), mm = amr_clock::mpi;
        if(depth>0 && depth<=MAXD)
        {
            const int k = stack[depth-1];
            t[k] += now-tlast;
            m[k] += mm-mlast;
        }
        tlast = now;
        mlast = mm;
    }
};

template<int N>
struct amr_phase
{
    amr_phase_clock<N> &c;
    bool on;
    amr_phase(amr_phase_clock<N> &cc, int k, bool root=false) : c(cc), on(root || cc.active()) { if(on) c.enter(k); }
    ~amr_phase() { if(on) c.leave(); }
    amr_phase(const amr_phase&) = delete;
    amr_phase& operator=(const amr_phase&) = delete;
};

#endif
