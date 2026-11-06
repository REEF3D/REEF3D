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

// Time step of the semi-implicit SFLOW (A 210 12/13, N 48 1): the advective Courant number of the flow
//   dt = N47 * 2 / max( |u|/dx + |v|/dy )            (same convention as sflow_etimestep; A 220 2/3: halved)
// plus the explicit diffusion limit (A 212 1), and the limit of the waves (cap):
//   A 228 > 0  the time step [s]
//   otherwise  from the grid, per wet cell with its width dx and water depth h, the smallest over the domain:
//              non-hydrostatic (A 220 2/3): a wave of 10 cells, lambda = 10 dx, with its phase speed c of the
//                model (omega^2 = g h k^2 (1 + B (kh)^2) / (1 + alpha/3 (kh)^2), B = (alpha - 1)/3, alpha A 224
//                for A 220 3, 1 for A 220 2) takes N steps per period: dt = lambda / (c N), N = 20 or -A 228;
//              hydrostatic (A 220 0): the gravity wave Courant number C: dt = C dx / sqrt(g h), C = 10 or -A 228.
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

public:
    // the wave limit without A 228 > 0 of a cell of width dx and water depth h (sflow_amr::si_timestep: the patches)
    static double wave_cap(lexer*, double dx, double h);

	double wd_criterion,cfl;
    double capw;        // the wave limit of the step (A 228 <= 0)
};

#endif
