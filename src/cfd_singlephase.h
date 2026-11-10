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

#ifndef CFD_SINGLEPHASE_H_
#define CFD_SINGLEPHASE_H_

#include"field1.h"
#include"field2.h"
#include"field3.h"
#include"field4.h"

class lexer;
class fdm;
class ghostcell;
class field;

// F 41 1: the single-phase level set. Only the water is in the momentum: the air takes the density
// and viscosity of the water (W 3 = W 1, W 4 = W 2, so that the face densities of the water at the
// surface are those of the water), the pressure is 0 at the free surface by the ghost fluid method
// of PIC single phase (pic_gfm.cpp: the Poisson rows, the gradients of the predictor and the
// correction, p = 0 in the air), and after every projection the faces between two cells in the air
// take the velocity of the water extended into them (3 layers, 0 beyond), as pic_f::air_velocity
// does for F 81 1. The level set and the convection near the surface then move with the water.
class cfd_singlephase
{
public:
    cfd_singlephase(lexer*);
    ~cfd_singlephase();

    void air_velocity(lexer*, fdm*, ghostcell*, field&, field&, field&);

    // the patches of the mesh refinement (cfd_amr::single_extend): the same extension component by
    // component, layer by layer over all patches of a level, the faces and the marks of the cells
    // around a patch from its siblings between the layers (the fills). ext_begin: the known faces
    // (the marks km), those around the patch known (the faces of the coarser grid there, extended
    // on it; the fill of the siblings sets their own marks); ext_layer: one layer; ext_end: 0 beyond
    void ext_begin(lexer*, fdm*, int);
    void ext_layer(lexer*, field&, int);
    void ext_end(lexer*, field&, int);
    field& marks(int c) { return c==0 ? static_cast<field&>(km1) : (c==1 ? static_cast<field&>(km2) : static_cast<field&>(km3)); }

    // F 41 2: the closed air pockets (cfd_pocket.h). pid: the pocket of a cell in the air at the last
    // labelling (>= 0), -1 a cell of the atmosphere, -2 none (water, a covered or a solid cell); pa: the
    // pressure of the pocket of a cell in the air, 0 in the atmosphere (the value of the ghost fluid
    // method at the surface, pic_gfm_air); lab: the labels of the labelling (scratch)
    field4 pid, pa, lab;

private:
    void extend(lexer*, fdm*, ghostcell*, field&, int, int);

    field1 km1;
    field2 km2;
    field3 km3;
    int fmarg;
};

#endif
