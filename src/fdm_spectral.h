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

#ifndef FDM_SPECTRAL_H_
#define FDM_SPECTRAL_H_

#include"increment.h"
#include"slice4.h"
#include"sliceint4.h"
#include"sliceint5.h"

class lexer;
class spectral_grid;
class spectral_store;

using namespace std;

/*--------------------------------------------------------------------
REEF3D::Spectral - field data on the 2D horizontal grid

  environment  bed, depth, eta, U, V: from the 2D grid in stand-alone
               runs, later from the coupled host model (SFLOW, NHFLOW)
  parameters   integrated wave parameters (spectral_param)
  wet          active cells: fluid (flagslice4 > 0) and depth >= A 705
  grid, N      spectral grid and block-sparse action density, owned by
               spectral_f
--------------------------------------------------------------------*/

class fdm_spectral : public increment
{
public:
    fdm_spectral(lexer*);
    virtual ~fdm_spectral() = default;

    // environment
    slice4 bed,depth,eta,U,V;

    // integrated wave parameters
    slice4 Hs,Tm01,Tm10,Tp,dir,spread;

    sliceint4 wet;
    sliceint5 nodeval;

    spectral_grid *grid = nullptr;
    spectral_store *N = nullptr;
};

#endif
