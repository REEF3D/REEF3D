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

#include"ship_models.h"
#include<cmath>

const double ship_models::Re_min = 1.0e5;

double ship_models::cf_ittc57(double Re)
{
    const double R = Re>Re_min ? Re : Re_min;
    const double lg = log10(R) - 2.0;
    
    return 0.075/(lg*lg);
}

double ship_models::friction(double rho, double nu, double S, double L, double k, double u, double &Re, double &CF)
{
    Re = fabs(u)*L/nu;
    CF = cf_ittc57(Re);
    
    return -0.5*rho*S*(1.0+k)*CF*fabs(u)*u;
}

void ship_models::crossflow(double rho, double Cd, const std::vector<double> &xs, const std::vector<double> &dx, const std::vector<double> &T,
                            double v, double r, double &Y, double &N)
{
    Y = N = 0.0;
    
    for(size_t s=0; s<xs.size(); ++s)
    {
        const double vl = v + xs[s]*r;   // local transverse velocity of the strip
        const double dY = -0.5*rho*Cd*T[s]*fabs(vl)*vl*dx[s];
        
        Y += dY;
        N += xs[s]*dY;
    }
}

double ship_models::roll_damping(double B44, double B44q, double p)
{
    return -B44*p - B44q*fabs(p)*p;
}
