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

//  CPM on the CFD mesh refinement (cfd_amr_cpm.cpp).
//
//  The parcels and the bed live on a grid of their own: the cells of the finest level (G 1) over the
//  region of the mobile bed (Q 85, else the bounding box of the boxes G 15, both widened to whole
//  level-0 cells), a regular grid of REEF3D decomposed over all ranks (its grid file built in memory,
//  decomposed by dmpart), with its own lexer, ghostcell and fdm, and the sediment objects of a CFD run
//  on it (sediment_part). CPM runs there unchanged. Every step the fluid of the finest grid at every
//  cell of it (velocities, pressure, level set, density, viscosity, eddy viscosity, k and its
//  dissipation) is copied onto it, from a patch of the finest level as it is, from a coarser grid
//  interpolated; its bed (topo) goes back onto every grid as the mean of the cells of CPM under each
//  cell.
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
    int px = 0, py = 0;         // periodic in x, y
    std::vector<int> pbox;      // the pieces of all ranks: origin and cells (6 per rank)

    // the fluid onto the grid of CPM: per rank the items requested from it (this rank: their position in
    // the arrays of the grid, sorted by rank) and served to it (the local grid, -1 level 0, and its cell:
    // the flat index of a cell of the finest level as it is, or a stencil of 8 cells per kind of field,
    // u v w and cell centred, on a coarser grid)
    std::vector<int> scnt, sdsp, rcnt, rdsp;
    std::vector<int> mine;
    std::vector<int> sv_g, sv_q, sv_st;         // sv_st: the stencil block, -1 a copy of sv_q
    std::vector<int> st_i;                      // stencil blocks: 4 kinds x 8 points (flat index, weight)
    std::vector<double> st_w;
    long ncopy = 0, nipol = 0;                  // cells of the grid copied, interpolated (this rank)
    long nipol_all = 0;                         //   interpolated (all ranks)
    std::vector<signed char> lev_of;            // per cell of this piece ((i knoy + j) knoz + k): the level of its fluid

    // the bed onto the grids: per cell of a grid in the region (a target, on
    // its rank) the cells of CPM under it on the ranks that hold them (items), summed there
    std::vector<int> bscnt, bsdsp, brcnt, brdsp;    // items: requested from / served to each rank
    std::vector<int> bt;                        // requested items: their target
    std::vector<int> tg, tq;                    // targets: grid (-1 level 0), flat index in its arrays
    std::vector<char> tleaf;                    //   an interior cell no finer patch covers (the bed is set)
    std::vector<double> tcnt;                   //   the cells of CPM under it
    std::vector<int> bsv;                       // served items: the box of cells of CPM (local, 6 per item)

    // G 27: the bed level of the columns at the start (this rank), the mobile columns of the grid of CPM
    // (all ranks, global column I*gn(L,1) + J of level L)
    std::vector<double> bed0;
    bool bed0_ok = false;
    std::vector<unsigned char> mob;
};

#endif
