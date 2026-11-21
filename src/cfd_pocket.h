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

#ifndef CFD_POCKET_H_
#define CFD_POCKET_H_

#include<vector>
#include<array>
#include<functional>

class lexer;
class fdm;
class ghostcell;
class momentum;
class cfd_amr;
class cfd_singlephase;

// F 41 2: the single-phase level set (F 41 1) with the closed air pockets. The air that is not
// connected to the atmosphere has a uniform pressure, adiabatic:
//
//   p = p_ref ((m / V)^gamma - 1)     (gauge; p_ref W 36, gamma W 33)
//
// with m the volume the air of the pocket would have at p_ref (its mass) and V its volume.
//
// Once per step, before the momentum, the cells in the air (phi < 0) of the leaf cells of all grids
// (not inside the solids and the topography) are labelled into connected components, across the
// ranks, the sibling patches and the patch edges.
// The component with the largest volume is the atmosphere (p = 0), as are those with fewer than 25
// (2D) or 125 (3D) cells in the air (not resolved).
//
// A pocket takes the mass and the volume of the pockets of the last labelling it overlaps (each
// shared by the overlaps of its cells), and the air of the atmosphere it overlaps at p_ref: a new
// pocket has p = 0, two pockets that merge add their masses, one that splits shares its mass.
// The volume of a pocket then changes with the flux of the water across its surface (the faces
// between its cells and the cells in the water, the velocity of the last step), drawn towards the
// volume of its level set by 2 % per step: the level set of a small pocket loses volume in the
// transport and the reinitialisation, which would read as compression; the flux alone drifts by
// some % within a second where the pocket is 8 cells across.
//
// The pressure follows its target from m and V with the factor 1/(1+G), G = gamma (p + p_ref) dt^2
// A^2 / (2 rho V^2), A the area of the surface in the water: the pressure of the pocket acts on the
// water of the next step, explicit. G estimates the gain of the loop pressure -> flux -> volume ->
// pressure over a step (the water around the pocket moved over a length 2 V / A); for G > 1 the
// explicit coupling would oscillate (small pockets, large steps), the factor makes it implicit
// there. A large pocket has G << 1. A pocket next to the surface of the atmosphere moves the water
// between them, a shorter length and a larger true gain: the coupling is not implicit then
// (validation 89, the bays of a deck).
//
// The volume of a component for the labelling: the cells in the air with the fraction 1/2 - phi/h
// (at most 1), and the cells in the water next to them with the same fraction, h the mean spacing.
//
// The pressure is the value at the surface in the ghost fluid method (pic_gfm_grad of the
// pressure, the predictor) and the pressure of the cells in the air (presscorr); the projection
// keeps it (its correction is 0 at the surface).
//
// The log: REEF3D_Log/REEF3D_CFD_Pockets.dat, the number of pockets and the volume of the level set,
// the volume and the pressure of the five largest per step.
// With the environment variable REEF3D_POCKET_DIAG set, rank 0 prints every pocket each step (its
// cells, volumes, flux, area, mass, pressures, gain and overlaps).

// a grid of the labelling: its cells in the list (i,j,k; the leaf cells, fluid)
struct cfd_pocket_grid
{
    lexer *p;
    fdm *a;
    cfd_singlephase *s;
    std::vector<std::array<int,3>> cells;
};

// two leaf cells next to each other across a patch edge: grid (index in the grid list) and cell;
// 1 the fine cell, 2 the coarse one at its side dir (+1, -1) along d; the velocity of the fine face
// (along +d) and its area
struct cfd_pocket_edge
{
    int g1, c1[3];
    int g2, c2[3];
    int d, dir;
    double uf, af;
};

// once per step before the momentum (driver_loop_cfd): level 0 alone or all grids of the mesh
// refinement (cfd_amr::pocket_grids)
void cfd_pocket_update(lexer*, fdm*, ghostcell*, momentum*, cfd_amr*);

// the labelling, the volumes and the pressures. exchange(k): the values of the ghost cells of all
// grids from their neighbours (the other ranks, the sibling patches; not across the patch edges),
// k 0 the labels (cfd_singlephase::lab), k 1 the pocket and its pressure (pid, pa)
void cfd_pocket_core(lexer*, std::vector<cfd_pocket_grid>&, const std::vector<cfd_pocket_edge>&, const std::function<void(int)>&);

#endif
