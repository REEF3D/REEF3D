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
class field;

// Particle-in-cell free surface flow for CFD (F 81; PICIN: Kelly et al. 2015, Chen et al. 2018)
// F 81 1: single-phase (water particles, p = 0 at the free surface)
// F 81 2: water particles, air phase Eulerian
// F 81 10: passive particles (phase 0)
//
// The water is seeded with particles (F 83 per cell and direction), which move with the CFD
// velocity, are exchanged between the subdomains and written to the LAGOON store (P 18).
// Phase 1 (F 81 1, 2): the free surface (phi) comes from the particles in place of the level set
// transport (pic_f_phi.cpp), particles are kept in the water (F 86, pic_f_reseed.cpp).
// Phase 2 (F 81 1, 2): the momentum of the water comes from the particles (APIC or FLIP, F 77,
// F 82; pic_f_flip.cpp). F 81 2: the air phase Eulerian, the pressure that of the two phases.
// F 81 1: single phase, no air in the momentum, p = 0 at the free surface (ghost fluid method,
// pic_gfm.cpp), the velocity of the water extended into the air (pic_f_extend.cpp).
// F 81 10: passive particles, the CFD unchanged.
//
// ini:   seeding in the water (phi > 0), phi from the particles, first output; after the
//        initial conditions, before the first CFD output
// step_begin (F 81 1, 2): the particle velocities to the faces in the water, before the momentum step
// start: every time step, after the momentum step and before the next time step is set:
//        advection with the new velocity over p->dt, walls and solids, MPI exchange,
//        phi from the particles and the density, keeping the particles in the water, output

class pic
{
public:
    virtual ~pic() = default;

    virtual void ini(lexer*, fdm*, ghostcell*)=0;
    virtual void start(lexer*, fdm*, ghostcell*)=0;

    // F 81 1, 2: the particle velocities to the water faces, before the momentum step
    virtual void step_begin(lexer*, fdm*, ghostcell*) {}

    // F 81 1, 2: the momentum step leaves out the convection on the faces in the water
    // (momentum_rk::convection_start; before: the right-hand side without the convection)
    virtual bool flip() const { return false; }
    virtual void mask_convection(lexer*, fdm*, int, const double*, field&) {}

    // F 81 1 (single phase, PICIN): no air in the momentum; after each projection the faces in the
    // air take the velocity of the water extended into them (momentum_rk::amr_project_after)
    virtual bool singlephase() const { return false; }
    virtual void air_velocity(lexer*, fdm*, ghostcell*, field&, field&, field&) {}
};

#endif
