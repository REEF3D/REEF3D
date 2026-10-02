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

#ifndef BC_NOFLUX_H_
#define BC_NOFLUX_H_

#include<vector>

class lexer;

/*--------------------------------------------------------------------
Faces of pressure-located cells (gcb4 list) where the normal flux through
the boundary is prescribed, so implicit operators must use a zero normal
gradient there (coefficient of the ghost cell into the diagonal) instead
of a lagged ghost value.

mask[IJK] has bit (cs-1) set for the face on side cs (1 x-, 2 y+, 3 y-,
4 x+, 5 z-, 6 z+).

  BC_NOFLUX_WALLS    lid (3), bed (5), walls and solid surfaces (21, 22)
  BC_NOFLUX_INFLOW   inflow / wave generation (1, 6) and patches with
                     prescribed velocity and Neumann pressure (2x1, 2x2 -> tens digit 1)
  BC_NOFLUX_SCALAR   leave out sides with a fixed scalar value (H 61-66),
                     for heat and concentration (ghost-cell label 61)
--------------------------------------------------------------------*/

enum
{
    BC_NOFLUX_WALLS  = 1,
    BC_NOFLUX_INFLOW = 2,
    BC_NOFLUX_SCALAR = 4
};

void bc_noflux_mask(lexer *p, std::vector<int> &mask, int what);

#endif
