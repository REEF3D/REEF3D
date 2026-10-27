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

#ifndef CFD_AMR_CPM_H_
#define CFD_AMR_CPM_H_

#include<vector>
#include<string>

class lexer;
class ghostcell;
class fdm;
class sediment;
class reinitopo;
class ioflow;
class solver;
class patchBC_interface;
class turbulence;
class field;

//  CPM on the CFD mesh refinement, static boxes (cfd_amr_cpm.cpp).
//
//  The parcels and the bed live on a grid of their own: the cells of the finest level (G 1) in the
//  bounding box of the refinement boxes G 15, a regular grid of REEF3D decomposed over all ranks
//  (its grid file built in memory, decomposed by dmpart), with its own lexer, ghostcell and fdm, and
//  the sediment objects of a CFD run on it (sediment_part). CPM runs there unchanged. Every step the
//  fluid of the finest patches (velocities, pressure, level set, density, viscosity, eddy viscosity,
//  k and its dissipation) is copied onto it, cell by cell, and its bed (topo) goes back onto the
//  patches, is restricted onto the coarser levels and sets their direct forcing.
struct cfd_amr_cpm
{
    lexer *p = nullptr;
    ghostcell *pgc = nullptr;
    fdm *a = nullptr;
    turbulence *pturb = nullptr;
    sediment *psed = nullptr;
    reinitopo *preto = nullptr;
    ioflow *pflow = nullptr;
    solver *psolv = nullptr;
    patchBC_interface *pBC = nullptr;

    int L = 0;                  // the level of the grid (the finest)
    int lo[3] = {0,0,0};        // its first cell on the index space of level L
    int n[3] = {0,0,0};         // its cells

    // the copy between the patches of level L and this grid: per rank the items requested from it
    // (this rank: their position in the arrays of the grid, sorted by rank) and served to it (the patch
    // and its cell); in the interior of the grid (the bed goes back from those only)
    std::vector<int> scnt, sdsp, rcnt, rdsp;    // items: sent to / received from each rank (requests)
    std::vector<int> mine;                      // requested items: the flat index in the arrays of the grid
    std::vector<char> inner;                    //   in the interior of the grid
    std::vector<int> sv_g, sv_q;                // served items: the local patch, the flat index in its arrays
    std::vector<char> sv_in;                    //   in the interior of the grid (the bed back)
};

#endif
