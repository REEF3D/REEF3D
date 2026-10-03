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

#ifndef DIFF_WALLGHOST_H_
#define DIFF_WALLGHOST_H_

#include"increment.h"
#include"field1.h"
#include<vector>

class lexer;
class ghostcell;

/*--------------------------------------------------------------------
Ghost-cell relation of a velocity component at the solid boundaries, for
an implicit treatment in the diffusion matrix.

At walls the ghost-cell routines set the first ghost cell as a linear
function of the interior cells on the wall-normal line,
    u(ghost) = c1 u(adjacent) + c2 u(next)
(no slip: Lagrange extrapolation through the wall value 0, slip: copy).
The implicit diffusion used the ghost value of the old stage field in the
rhs, so the wall flux lagged one stage behind. With c1, c2 the ghost cell
is coupled to the new values instead (c1 into the diagonal, c2 into the
coefficient of the next cell).

c1 and c2 are obtained from the ghost-cell routines themselves (probe
fields with 0, 1, m and m^2, m = i+j+k global), so they follow whatever
B 20 / B 23 / B 29 select. A face is only treated implicitly if its ghost
value is homogeneous (no prescribed value, e.g. inflow), depends only on
those two cells, and the ghost cell belongs to this face alone.
--------------------------------------------------------------------*/

class diff_wallghost : public increment
{
public:
    diff_wallghost(lexer*);
    virtual ~diff_wallghost();

    // c = 0, 1, 2 (u, v, w), gcv: ghost-cell label of the velocity component
    void update(lexer*, ghostcell*, int c, int gcv);

    // face cs (1 x-, 2 y+, 3 y-, 4 x+, 5 z-, 6 z+) of cell (ii,jj,kk), after update() for its component
    bool coef(lexer*, int ii, int jj, int kk, int cs, double &c1, double &c2);

private:
    void probe(lexer*, ghostcell*, int, int*, int, std::vector<double>&);
    int ghost_index(lexer*, int, int, int, int);

    field1 P;
    std::vector<int> owner;
    std::vector<double> w1, w2;
    std::vector<char> valid;
    int **gcb;
    int gcb_count;
    int gcv_probe;
};

#endif
