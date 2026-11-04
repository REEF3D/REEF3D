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

#ifndef PIC_H_
#define PIC_H_

class lexer;
class fdm;
class ghostcell;

// Particle-in-cell free surface flow for CFD (F 81; PICIN: Kelly et al. 2015, Chen et al. 2018)
// F 81 1: single-phase (water particles, p = 0 at the free surface)
// F 81 2: water particles, air phase Eulerian
//
// Phase 0: the particle infrastructure only. The water is seeded with particles (F 83 per
// cell and direction), which are moved passively with the CFD velocity, exchanged between
// the subdomains and written to the LAGOON store (P 18). The CFD solution is not changed.
//
// ini:   seeding in the water (phi > 0), first output; after the initial conditions
// start: every time step, after the momentum step and before the next time step is set:
//        advection with the new velocity over p->dt, walls and solids, MPI exchange, output

class pic
{
public:
    virtual ~pic() = default;

    virtual void ini(lexer*, fdm*, ghostcell*)=0;
    virtual void start(lexer*, fdm*, ghostcell*)=0;
};

#endif
