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


#ifndef SEASTATE_NHFLOW_H_
#define SEASTATE_NHFLOW_H_

#include"seastate_coupling.h"

class fdm_nhf;
class nhflow_amr;

using namespace std;

/*--------------------------------------------------------------------
REEF3D::SEASTATE <-> NHFLOW coupling (A 10 5, A 750 1)

The NHFLOW host of seastate_coupling (formulations there). SEASTATE runs
on the 2D slice grid of NHFLOW; the layered host sees the waves through
depth-averaged quantities:

  start        every NHFLOW step, before the momentum step: depth-
               averaged velocity Ub = sum U DZN, environment (eta, WL,
               Ub, Vb), wave step when due, wave force
  u_source     wave force [m^2/s^2] in the momentum right-hand side of
  v_source     every layer (uniform over the depth: sum F DZN = F_w),
               every RK stage
  mass_source  -div(M) in the continuity of every RK stage (vortex
               force, A 751 2)
  ghostcells   surfbeat (A 770 1): long-wave ghost cells at the x- side,
               depth-uniform u_g in every layer, v and w of the first
               cell (zero gradient); called after the inflow boundary
               conditions of the step for the NHFLOW and the RK arrays
  wl_ghostcells  WL and eta of the ghost cells: zero gradient
               (Flather: h_g = h_1)

With surfbeat the x- edge is open for the NHFLOW ghost cells (lexer
open_xm = 1, as the Flather/Riemann edges of iowave), so that ghostcell
keeps the long-wave values. Not with mesh refinement (G 1).

Curvilinear and multi-block grids (G 50): the remote waves only (A 750 2,
SEASTATE on its own grid). The host maps between its grids and the wave
grid at the physical cell centres (seastate_coupling_grids.cpp); on a
multi-block grid every block has a seastate_nhflow of its own that holds
its wave fields (Fw, Gw, divM, Hs, Tp; no wave model), so the momentum of
the block applies them as level 0 does. The set-up exchange waits for the
blocks (attach, after nhflow_amr::ini).
--------------------------------------------------------------------*/

class seastate_nhflow : public seastate_coupling
{
public:
    seastate_nhflow(lexer*, fdm_nhf*, ghostcell*);
    virtual ~seastate_nhflow() = default;

    void ini(lexer*, fdm_nhf*, ghostcell*);

    // multi-block grid (G 50 1): the blocks of nhflow_amr (after its ini) get their seastate_nhflow, then the set-up
    // exchange with the remote waves
    void attach(lexer*, fdm_nhf*, ghostcell*, nhflow_amr*);

    // every NHFLOW step, before the momentum step: wave step when due
    void start(lexer*, fdm_nhf*, ghostcell*);

    // momentum right-hand side (every RK stage, every layer) and continuity source (A 751 2)
    void u_source(lexer*, fdm_nhf*);
    void v_source(lexer*, fdm_nhf*);
    void mass_source(lexer*, fdm_nhf*, slice&);

    // long-wave boundary (surfbeat): ghost cells of the x- side
    void ghostcells(lexer*, fdm_nhf*, double *UH, double *VH, double *WH);
    void wl_ghostcells(lexer*, fdm_nhf*, slice &WL);

private:
    void average(lexer*, fdm_nhf*);

    // the grids of a curvilinear or multi-block host on this rank (level 0 first), their depth-averaged velocities
    void host_grids(lexer*, fdm_nhf*, std::vector<seastate_host_grid>&);

    slice4 Ub,Vb;

    nhflow_amr *pamr = nullptr;     // multi-block grid: the blocks
    bool pending = false;           // the set-up exchange waits for the blocks (attach)
};

#endif
