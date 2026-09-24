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

#ifndef SFLOW_HLL_H_
#define SFLOW_HLL_H_

#include"increment.h"

class lexer;
class fdm2D;
class slice;
class ghostcell;
class patchBC_interface;
class sflow_flux_build;

using namespace std;

// HLL finite volume fluxes for the depth-averaged SFLOW equations on the
// cell-centred slice grid (conserved variables WL, UH, VH, WH).
//   ipol 1: UH  -> b->F     ipol 2: VH -> b->G
//   ipol 3: WH  -> b->H     ipol 4: continuity -> b->FEx, b->FEy
// The momentum fluxes contain the hydrostatic part 0.5*g*eta^2 + g*eta*d_face,
// the matching bed-slope source is added by sflow_pressure::upgrad/vpgrad.

class sflow_HLL final : public increment
{
public:
	sflow_HLL(lexer*,ghostcell*,patchBC_interface*);
	virtual ~sflow_HLL();

    void start(lexer*, fdm2D*, int);

private:
    void aij_U(lexer*, fdm2D*);
    void aij_V(lexer*, fdm2D*);
    void aij_W(lexer*, fdm2D*);
    void aij_E(lexer*, fdm2D*);

    void HLL(lexer*, fdm2D*, slice&, slice&, slice&, slice&, slice&, slice&);
    void flux_bc(lexer*, fdm2D*, int);
    void divergence(lexer*, fdm2D*, slice&);

    double denom;
    int inflow, outflow;

    ghostcell *pgc;
    patchBC_interface *pBC;
    sflow_flux_build *pflux;
};

#endif
