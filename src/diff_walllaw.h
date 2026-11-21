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

#ifndef DIFF_WALLLAW_H_
#define DIFF_WALLLAW_H_

#include"increment.h"
#include<unordered_set>

class lexer;
class matrix_diag;

/*--------------------------------------------------------------------
The faces of a velocity component where the wall law of the momentum equations (B 10, bcmom) puts
the wall stress: the walls of the ghost-cell lists gcb of the types 5, 21, 22, 41, 42, 43 (not on a
periodic boundary) and the direct-forcing faces gcdf (solid, topography, 6DOF bodies), tangential
to the component. The same selection as bcmom::bcmom_start.

With B 10 != 0 the diffusion has zero gradient through these faces, so that the wall law is the only
wall stress there (D 22 1 does that for all faces to walls and direct-forcing cells; D 22 2, the PLIC
diffusion and the explicit diffusion took the gradient to the ghost or solid cell as well).
--------------------------------------------------------------------*/

class diff_walllaw : public increment
{
public:
    diff_walllaw();

    // c = 0, 1, 2 (u, v, w); rebuilds the face set from the current lists
    void update(lexer*, int c);

    // face cs (1 x-, 2 y+, 3 y-, 4 x+, 5 z-, 6 z+) of cell (ii,jj,kk) carries the wall law
    bool face(lexer*, int ii, int jj, int kk, int cs) const;

    // zero gradient through these faces in the rows of the implicit diffusion matrix (D 20 2)
    void rows(lexer*, matrix_diag&, int c);

    bool active;

private:
    long long key(lexer*, int, int, int, int) const;
    std::unordered_set<long long> S;
};

#endif
