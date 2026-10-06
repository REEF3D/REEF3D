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


#ifndef SEASTATE_SFLOW_H_
#define SEASTATE_SFLOW_H_

#include"seastate_coupling.h"

class fdm2D;

using namespace std;

/*--------------------------------------------------------------------
REEF3D::SEASTATE <-> SFLOW coupling (A 10 2, A 750 1)

The SFLOW host of seastate_coupling (formulations there):

  start        every SFLOW step, before the momentum step: environment
               (SFLOW eta, WL, U, V), wave step when due, wave force
  u_source     wave force in the SFLOW momentum right-hand side, every
  v_source     RK stage
  mass_source  -div(M) in the continuity (vortex force, A 751 2)
  ghostcells   surfbeat (A 770 1): long-wave ghost cells at the x- side
               (called in sflow_momentum_func::ghostcells after the
               velocity boundary conditions, so that WL, UH, VH of the
               ghost cells follow)
  flux_bc      surfbeat: fluxes through the x- faces from the ghost
               state (called at the end of sflow_HLL::flux_bc)
--------------------------------------------------------------------*/

class seastate_sflow : public seastate_coupling
{
public:
    seastate_sflow(lexer*, fdm2D*, ghostcell*);
    virtual ~seastate_sflow() = default;

    void ini(lexer*, fdm2D*, ghostcell*);

    // every SFLOW step, before the momentum step: wave step when due
    void start(lexer*, fdm2D*, ghostcell*);

    // SFLOW momentum right-hand side (every RK stage) and continuity source (A 751 2)
    void u_source(lexer*, fdm2D*);
    void v_source(lexer*, fdm2D*);
    void mass_source(lexer*, fdm2D*, slice&);

    // long-wave boundary (surfbeat): ghost cells of the x- side, and the x- face fluxes
    // (physical flux of the ghost state, as at SFLOW inflow boundaries; SFLOW takes the
    // other boundary faces as walls)
    void ghostcells(lexer*, fdm2D*);
    void flux_bc(lexer*, fdm2D*, int ipol, slice &Fx);
};

#endif
