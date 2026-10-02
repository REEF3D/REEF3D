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

#ifndef SFLOW_AMR_LINESOLVE_H_
#define SFLOW_AMR_LINESOLVE_H_

#include"solver2D.h"
#include"increment.h"
#include<vector>

using namespace std;

//  Solver for the line-implicit systems of a refined SFLOW patch (the u_a inversion of the
//  Boussinesq equations, A 220 4): the rows couple along x only or along y only and are solved
//  exactly line by line (Thomas); the values beyond the computed cells of the patch are the
//  filled ghost values of f (Dirichlet).  Other matrices: Gauss-Seidel sweeps.

class sflow_amr_linesolve final : public solver2D, public increment
{
public:
    void start(lexer*, ghostcell*, slice&, matrix2D&, vec2D&, vec2D&, int) override final;

private:
    vector<int> id;
    vector<double> a,bb,c,d;
};

#endif
