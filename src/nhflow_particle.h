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

#ifndef NHFLOW_PARTICLE_H_
#define NHFLOW_PARTICLE_H_

class lexer;
class fdm_nhf;
class ghostcell;

// Lagrangian particle tracking for NHFLOW (passive tracers, floating and buoyant macroplastics)
// step_begin: release + velocity at (x^n,t^n), before the momentum step
// step_end:   Heun corrector with the field at t^n+1, physics, MPI exchange, output

class nhflow_particle
{
public:
    virtual ~nhflow_particle() = default;

    virtual void ini(lexer*, fdm_nhf*, ghostcell*)=0;
    virtual void step_begin(lexer*, fdm_nhf*, ghostcell*)=0;
    virtual void step_end(lexer*, fdm_nhf*, ghostcell*)=0;
};

#endif
