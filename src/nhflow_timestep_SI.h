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

#ifndef NHFLOW_TIMESTEP_SI_H_
#define NHFLOW_TIMESTEP_SI_H_

#include"nhflow_timestep.h"

using namespace std;

// time step of the hydrostatic semi-implicit NHFLOW (A 534 1): the advective Courant number
//   dt = 2 N47 / max( |u|/dx + |v|/dy + |omega|/(dsig WL) ),
// the limit of the explicit horizontal diffusion, capped by A 536 (s); without A 536 the time step stays
// at the explicit gravity wave limit
class nhflow_timestep_SI final : public nhflow_timestep
{
public:
	nhflow_timestep_SI(lexer*);
	virtual ~nhflow_timestep_SI();
	void start(lexer*,fdm_nhf*,ghostcell*) override;
	void ini(lexer*,fdm_nhf*,ghostcell*) override;

private:
    const double cfl;
};

#endif
