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
// REEF3D::CORAL: steady operating curve of a rotor and the data for its controller, from CORAL's own blade element
// momentum (rotor_steady: axial inflow, no cone or tilt), as the ROSCO toolbox does with its Cp surface:
//   Cp(lambda) at zero pitch: Cp_max at lambda_opt, the region-2 torque gain K = rho pi R^5 Cp_max/(2 lambda_opt^3)
//   (rotor side, Q = K Omega^2; the generator side K/N^3)
//   the operating curve over the wind speed: below rated Omega = lambda_opt U/R up to Omega_rated at zero pitch; the
//   rated wind speed where the power reaches P_rated; above it the pitch that holds P_rated at Omega_rated and the power
//   sensitivity dP/dtheta there, fitted as dP/dtheta(theta) = dP/dtheta(0) (1 + theta/theta_K) (gain scheduling) over
//   the points above 1 deg of pitch; if the fitted line is not negative at zero pitch (the power rises when the rotor
//   starts to pitch at rated), dP/dtheta(0) is the line's value at the first fitted pitch (lower, safe gains)

#ifndef CORAL_ROTOR_TUNE_H_
#define CORAL_ROTOR_TUNE_H_

#include"coral/aero/coral_bem.h"
#include<string>
#include<vector>

namespace coral
{
    struct rotor_tuning
    {
        double cp_max = 0.0, lambda_opt = 0.0, K = 0.0;      // K rotor side [N m s^2]
        double U_rated = 0.0;
        double dPdtheta0 = 0.0, theta_K = 0.0;              // W/rad, rad
        struct point {double U, Omega, theta, P, T, Cp, dPdtheta;};
        std::vector<point> curve;                          // 3 m/s .. U_out

        bool write(const std::string &file) const;
    };

    // rotor at rated power P_rated (aerodynamic) and rated rotor speed Omega_rated [rad/s]
    rotor_tuning tune_rotor(const rotor_def&, double rho, double P_rated, double Omega_rated, double U_out=25.0,
                            const bem_options& = bem_options());
}

#endif
