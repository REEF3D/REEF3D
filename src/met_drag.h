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
       5  Garratt capped at W 213 (default 2.5e-3): the saturation of
          the drag in hurricane winds (Powell et al. 2003, Donelan et
          al. 2004), reached at 26.1 m/s; the default of the vortex W 200 3
       6  Zijlema et al. (2012), beyond its peak (2.03e-3 at 31.4 m/s)
          not below W 214 (default 1.5e-3): the decline of Powell et
          al. (2003) without the unphysical zero at 68 m/s
W 207 0 (default): 1 for W 200 1/2 (as before), 5 for the vortex and the
          blend (W 200 3, 4).

Laws 1-4 are those of wind_f (A 570), evaluated per cell with the local
wind speed.
--------------------------------------------------------------------*/

// the law of W 207 (0: auto)
inline int met_drag_law(int W207, int W200) {return W207>0 ? W207 : (W200>=3 ? 5 : 1);}

inline double met_drag(int law, double U, double cap=2.5e-3, double floor=1.5e-3)
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

    if(law==5)
    {
    Cd = 0.001*(0.75 + 0.067*U);
    Cd = Cd<cap ? Cd : cap;
    }

    if(law==6)
    {
    const double x = U/31.5;
    Cd = 0.001*(0.55 + 2.97*x - 1.49*x*x);

        if(x>2.97/2.98)
        Cd = Cd>floor ? Cd : floor;
    }

    return Cd;
}

#endif
