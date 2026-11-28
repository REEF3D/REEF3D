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

#ifndef PLIC_DISTANCE_H_
#define PLIC_DISTANCE_H_

#include"field4.h"

class lexer;
class ghostcell;

// The signed distance (water > 0) to the PLIC interface of a VOF field, for the output of a
// VOF-PLIC run (F 80 4) only: the printed phi. The solver's own phi of a PLIC run is a phase
// marker (+1 water, -1 air, (vof - 0.5) dz in the interface cells, VOF_PLIC::start), whose zero
// lies at the cell centres; this one passes through the PLIC planes.
//
//  - every interface cell (F 93 <= vof <= F 94): Youngs normal from the 3x3(x3) vof stencil and
//    the plane constant of Scardovelli & Zaleski (as VOF_PLIC::reconstructPlane_alt); a full cell
//    next to an empty one across a face: the plane of that face
//  - the cells within two cells of such a plane: the shortest distance to the facets near them
//    (the plane clipped to its cell),
//    with the sign of their own vof - 0.5
//  - further away: the vertical distance to the free surface of the column (top face of the
//    highest fluid cell minus the air below it), at least the band, with the sign of vof - 0.5
//
// Local to the rank (the planes of the ghost cells too), then one ghost cell exchange. Computed
// at the output times only; the solution is not changed.
class plic_distance
{
public:
    plic_distance(lexer*);

    field4& start(lexer*, ghostcell*, field &vof);

    // the plane constant of a cell (dx, dy, dz) with the fraction V0 below the plane of unit
    // normal (nx, ny, nz), pointing into the air: n.(x - centre) = alpha
    static double plane_constant(double nx, double ny, double nz, double dx, double dy, double dz, double V0);

private:
    void candidate(lexer*, int ii, int jj, int kk, double nx, double ny, double nz, double px, double py, double pz);

    field4 dist;
    const int band = 2;
};

#endif
