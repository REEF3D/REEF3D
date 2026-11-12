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

#ifndef SI_TRIDIAG2D_H_
#define SI_TRIDIAG2D_H_

#include"solver2D.h"
#include"increment.h"
#include<vector>
#include<memory>

class lexer;
class ghostcell;
class slice;
class matrix2D;
class vec2D;

// The eta system of the semi-implicit NHFLOW with one row of cells (j_dir 0) on one rank: tridiagonal in the
// SLICELOOP4 order, solved directly (Thomas). Rows coupled to a non-neighbouring row (or a coupling in y): the
// iterative solver of the slice (fallback) for that solve.
class si_tridiag2D : public solver2D, public increment
{
public:
    si_tridiag2D(lexer*, ghostcell*, solver2D *fallback);
    virtual ~si_tridiag2D();

    void start(lexer*, ghostcell*, slice&, matrix2D&, vec2D&, vec2D&, int) override;

private:
    solver2D *fb;
    std::vector<double> c, d;
    std::vector<int> ci;
};

#endif
