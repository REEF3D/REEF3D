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

#ifndef MATRIX_DIAG_H_
#define MATRIX_DIAG_H_

#include <vector>

class lexer;

class matrix_diag
{
public:
    matrix_diag(lexer*);

    //  Explicit row count, for a solver path that knows it writes fewer rows
    //  than lexer's veclength allows for - see fdm_fnpf.
    matrix_diag(lexer*, int rows);
    virtual ~matrix_diag() = default;

    void resize(int);

    void reset();

    std::vector<double> n,s,e,w,b,t,p;

    //  Optional x-sigma and y-sigma couplings of the FNPF Laplace operator
    //  (A 328 1): nt = (i+1,k+1), nb = (i+1,k-1), st = (i-1,k+1), ...,
    //  w = j+1, e = j-1. Empty unless the Laplace assembly sizes them; only
    //  REEFMG's Krylov operator reads them.
    std::vector<double> nt,nb,st,sb,wt,wb,et,eb;
};

#endif
