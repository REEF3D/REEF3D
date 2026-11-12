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

#ifndef CFD_GRAVITYWAVE_DT_H_
#define CFD_GRAVITYWAVE_DT_H_

#include<cmath>

// N 54 > 0: the gravity-wave limit of the CFD time step. Without gravity in the acceleration term the
// step must still resolve the shortest gravity wave of the grid at the free surface: wavelength 2 dx,
// k = pi/dx, deep water omega = sqrt(g k). dt <= N54/omega; 1.0e10 when N 54 is off or there is no gravity.
inline double cfd_gravitywave_dt(double N54, double gx, double gy, double gz, double dxmin)
{
    const double g = std::sqrt(gx*gx + gy*gy + gz*gz);

    if(N54<=0.0 || g<=0.0 || dxmin<=0.0 || dxmin>=1.0e9)
    return 1.0e10;

    return N54/std::sqrt(M_PI*g/dxmin);
}

#endif
