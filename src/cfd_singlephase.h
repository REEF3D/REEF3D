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

private:
    void extend(lexer*, fdm*, ghostcell*, field&, int, int);

    field1 km1;
    field2 km2;
    field3 km3;
    int fmarg;
};

#endif
