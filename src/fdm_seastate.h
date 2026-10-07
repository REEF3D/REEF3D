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

#ifndef FDM_SEASTATE_H_
#define FDM_SEASTATE_H_

#include"increment.h"
#include"slice4.h"
#include"sliceint4.h"
#include"sliceint5.h"

class lexer;
class seastate_grid;
class seastate_store;

using namespace std;

/*--------------------------------------------------------------------
REEF3D::SEASTATE - field data on the 2D horizontal grid

  environment  bed, depth, eta, U, V: from the 2D grid in stand-alone
               runs, later from the coupled host model (SFLOW, NHFLOW)
  parameters   integrated wave parameters (seastate_param)
  wet          active cells: fluid (flagslice4 > 0), inside the global
               domain and depth >= A 705
  wet0         the active cells at the start, which have spectral storage;
               a coupled host can dry and re-wet them (wet <= wet0)
  kinematics   depth and current gradients, depth of the previous step;
               refr = 1 where no neighbour is land or dry (refraction and
               frequency shift are switched off next to dry cells, as in
               SWAN; at the edge of the domain one-sided gradients)
  grid, N      spectral grid and block-sparse action density
  kw, cg       wave number and group velocity per active cell and
               frequency (block-sparse, nbin = nsig)
  all owned by seastate_f
--------------------------------------------------------------------*/

class fdm_seastate : public increment
{
public:
    fdm_seastate(lexer*);
    virtual ~fdm_seastate() = default;

    // environment
    slice4 bed,depth,eta,U,V;

    // integrated wave parameters
    slice4 Hs,Tm01,Tm10,Tp,dir,spread;

    // kinematics
    slice4 ddx,ddy,dUdx,dUdy,dVdx,dVdy,dddt,depth_n;

    sliceint4 wet,wet0,refr;
    sliceint5 nodeval;

    seastate_grid *grid = nullptr;
    seastate_store *N = nullptr;
    seastate_store *kw = nullptr;
    seastate_store *cg = nullptr;
};

// k, cg, depth and current gradients, refraction flags (seastate_f_kinematics.cpp); (i+gi0, j+gj0)
// is the global index of local cell (i,j), the domain has gnx x gny cells of the grid's level
void seastate_kinematics(lexer*, fdm_seastate*, double dt_depth, int gi0, int gj0, int gnx, int gny);

#endif
