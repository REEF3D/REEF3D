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

#ifndef NHFLOW_SCOUR_H_
#define NHFLOW_SCOUR_H_

#include"scour_core.h"

class fdm_nhf;

// P 187, P 198: the scour parameters of NHFLOW (scour_core.h). The near-bed velocity is the velocity of
// the lowest sigma layer; the computed current is the time mean of the depth-averaged velocity, added
// to the ambient current of P 189 or P 196 (default zero). U_m is the orbital part of the near-bed
// velocity (its standard deviation), so the mean current near the bed does not enter it.
class nhflow_scour : public scour_core
{
public:
    nhflow_scour(lexer*,fdm_nhf*,ghostcell*);
    virtual ~nhflow_scour() = default;

    // every time step: the near-bed and depth-averaged velocity and their statistics
    void sample(lexer*,fdm_nhf*,ghostcell*);
};

#endif
