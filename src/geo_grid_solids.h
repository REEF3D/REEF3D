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

#ifndef GEO_GRID_SOLIDS_H_
#define GEO_GRID_SOLIDS_H_

#include<vector>

class lexer;
struct geo_cart;

// the inside (-1) / outside (1) flags of the entities of a role (geo_mesh::role_solid, role_topo)
// of the grid file on the cells of g of lexer p, as lexer::grid_solids casts them (cells outside the
// global domain: outside); false: no entities of the role.  The CFD mesh refinement casts the solids
// and the topography on its patches with it (cfd_amr_geom.cpp)
bool geo_grid_solids_io(lexer *p, int role, const geo_cart &g, std::vector<int> &IO);

#endif
