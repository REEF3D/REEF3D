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

#ifndef SFLOW_EDIFF_H_
#define SFLOW_EDIFF_H_

#include"sflow_diffusion.h"
#include"increment.h"

class lexer;
class fdm2D;
class ghostcell;
class solver2D;
class slice;
class sliceint;

using namespace std;

class sflow_ediff final : public sflow_diffusion, public increment
{
public:
	sflow_ediff(lexer*);

    // the grid computes its eddy viscosity itself (A 260 3 / 4) when diff_u is called: curvilinear grids and the
    // grids of REEFAMR, whose state the turbulence module of the native grid does not see (sflow_amr, sflow_f_logic)
    int own_eddyv = 0;
	virtual ~sflow_ediff();

	void diff_u(lexer*, fdm2D*, ghostcell*, solver2D*, slice&, slice&, slice&, slice&, slice&, double) override final;
	void diff_v(lexer*, fdm2D*, ghostcell*, solver2D*, slice&, slice&, slice&, slice&, slice&, double) override final;
    void diff_w(lexer*, fdm2D*, ghostcell*, solver2D*, slice&, slice&, slice&, slice&, slice&, slice&, double) override final;
    void diff_scalar(lexer*, fdm2D*, ghostcell*, solver2D*, slice&, double, double) override final;

    // the algebraic eddy viscosity (A 260 3 Prandtl mixing length, A 260 4 parabolic) of the cells of the grid
    // and of its first ghost layer
    static void algebraic_eddyv(lexer*, fdm2D*, slice &U, slice &V, slice &WL);

private:
    double viscosity(lexer*, fdm2D*);
    // curvilinear grids: the depth-integrated stress through the faces of the cells, comp 1 x (to F), 2 y (to G)
    void curv_stress(lexer*, fdm2D*, slice &U, slice &V, slice &WL, int comp);
    double dxx(lexer*, slice&);
    double dyy(lexer*, slice&);
    double dxy(lexer*, slice&);
};

#endif
