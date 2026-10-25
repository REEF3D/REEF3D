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

#ifndef MET_DRAG_H_
#define MET_DRAG_H_

#include<cmath>

/*--------------------------------------------------------------------
REEF3D atmospheric forcing (W 200): drag coefficient of the 10 m wind,
tau = rho_air Cd(|U10|) |U10| U10.  No lexer, unit test met_test.cpp.

W 207  1  Garratt (1977), capped at 3.5e-3 above 41 m/s (as A 570 1)
       2  Wu (1982)                                        (as A 570 2)
       3  Smith and Banke (1975)                           (as A 570 3)
       4  Zijlema et al. (2012), Uref = 31.5 m/s, >= 0     (as A 570 4)

The laws are those of wind_f (A 570), evaluated per cell with the local
wind speed. Saturating high-wind laws for tropical cyclones follow in a
later patch.
--------------------------------------------------------------------*/

inline double met_drag(int law, double U)
{
    double Cd = 0.0;

    if(law==1)
    Cd = (U<41.0) ? 0.001*(0.75 + 0.067*U) : 0.001*3.5;

    if(law==2)
    Cd = (U<7.5) ? 1.2875e-3 : 0.001*(0.8 + 0.065*U);

    if(law==3)
    Cd = 0.001*(0.63 + 0.066*U);

    if(law==4)
    {
    const double x = U/31.5;
    Cd = 0.001*(0.55 + 2.97*x - 1.49*x*x);
    Cd = Cd>0.0 ? Cd : 0.0;
    }

    return Cd;
}

#endif
