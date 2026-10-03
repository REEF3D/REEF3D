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

#include"CPM.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment_fdm.h"
#include"turbulence.h"

// Andrews & O'Rourke (1996):  Dp = Cd 3/8 rho_f/rho_p |Uf-Up|/r_p,  Cd = 24/Re (thf^-2.65 + 1/6 Re^2/3 thf^-1.78)
// written without the 1/Re singularity:  Dp = 18 nu rho_f/(rho_p d^2) (thf^-2.65 + 1/6 Re^2/3 thf^-1.78)
// vel: magnitude of the relative velocity |Uf-Up|
double CPM::drag_model(lexer *p, double d, double rhoS, double vel, double Ts)
{
    double Tf = MAX(1.0-Ts, 1.0-theta_max);
    Tf = MIN(Tf,1.0);

    double Rep = fabs(vel)*d/p->W2;

    double Dp = 18.0*p->W2*p->W1/(rhoS*d*d) * (pow(Tf,-2.65) + (1.0/6.0)*pow(Rep,2.0/3.0)*pow(Tf,-1.78));

    return Dp;
}
