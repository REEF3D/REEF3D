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

#ifndef SFLOW_TIMESTEP_SI_H_
#define SFLOW_TIMESTEP_SI_H_

#include"sflow_timestep.h"
#include"increment.h"

using namespace std;

// Time step of the semi-implicit SFLOW (A 226 1, N 48 1): the advective Courant number
//   dt = N47 * 2 / max( |u|/dx + |v|/dy )            (same convention as sflow_etimestep)
// plus the explicit diffusion limit (A 212 1), capped by A 228 (maximum time step, s).
// Without A 228 the cap is the explicit gravity wave limit (no gain, but safe).
// The gravity wave Courant number dt sqrt(g h_max)/dx_min is printed.

class sflow_timestep_SI final : public sflow_timestep, public increment
{
public:
	sflow_timestep_SI(lexer*,fdm2D*);
	virtual ~sflow_timestep_SI();

    void start(lexer*,fdm2D*,ghostcell*) override final;
	void ini(lexer*,fdm2D*,ghostcell*) override final;

private:
    double advective(lexer*,fdm2D*,ghostcell*,double&);

	double wd_criterion,cfl;
};

#endif
