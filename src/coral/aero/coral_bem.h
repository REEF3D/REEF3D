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
// REEF3D::CORAL: blade element momentum theory.
//
// One blade element (bem_solve): the inflow at the element without induction, Vx normal to the rotor plane
// (downwind) and Vy in the plane (the blade's rotation, positive in the direction it moves), give the inflow angle
// phi with
//   tan(phi) = Vx (1 - a) / (Vy (1 + a')),   alpha = phi - theta   (theta: twist + pitch, towards feather)
//   cn = cl cos(phi) + cd sin(phi),  ct = cl sin(phi) - cd cos(phi),  sigma' = B c / (2 pi r)
//   k = sigma' cn / (4 F sin^2 phi),  k' = sigma' ct / (4 F sin phi cos phi)
//   a = k/(1 + k) for k <= 2/3, else Buhl's empirical high-induction relation (C_T = 8/9 + (4F - 40/9) a + (50/9 - 4F) a^2),
//   a' = k'/(1 - k'); propeller brake (phi < 0): a = k/(k - 1) for k > 1, else 0
//   F = F_tip F_hub, Prandtl: F_tip = 2/pi acos(exp(-B (R - r) / (2 r |sin phi|))), F_hub likewise with (r - R_hub)/R_hub
// solved for phi alone with Ning's residual (Wind Energy 17, 2014; as AeroDyn 15): bracketed in (0, pi/2], then the
// propeller brake (-pi/4, 0), then (pi/2, pi), Brent's method; where the ends of none of them change sign (a pole of
// a' inside, e.g. at the root of a blade with a large chord), the three intervals are scanned in 200 parts for a
// sign change that is a root; without one: no induction and converged = false.
// Loads per unit span: fn = 1/2 rho W^2 c cn (along Vx, thrust), ft = 1/2 rho W^2 c ct (along Vy, driving torque),
// mz = 1/2 rho W^2 c^2 cm (pitching moment, nose up positive), W^2 = (Vx (1 - a))^2 + (Vy (1 + a'))^2.
// Without inflow normal to the plane or in it (Vx <= 0 or Vy <= 0), or with induction off: no induction. At the tip
// (r >= R) with the tip loss, and at the hub (r <= R_hub) with the hub loss, F = 0: no load.
// Skewed wake (Pitt-Peters, AeroDyn's SkewMod 1): a (1 + skew) after the solution, with
// skew = 15 pi/32 (r/R) tan(chi/2) cos(psi) given by the caller (chi the wake skew, psi the azimuth from the most
// downwind position of the blade).
//
// rotor_steady: the rotor in uniform axial inflow U at the speed Omega and pitch, rigid and straight blades (no cone,
// tilt or yaw): the elements at the blade stations, thrust and torque by the trapezoid rule over the stations.

#ifndef CORAL_BEM_H_
#define CORAL_BEM_H_

#include"coral/aero/coral_airfoil.h"
#include<string>
#include<vector>

namespace coral
{
    struct bem_options
    {
        bool induction = true;
        bool tip_loss = true, hub_loss = true;
        bool tangential = true;        // a'
        bool drag_in_induction = true; // cd in cn and ct of the momentum balance (AeroDyn's AIDrag, TIDrag)
        bool buhl = true;              // Buhl's high-induction correction, else the momentum relation a = k/(1 + k) for all k
        double tol = 1.0e-10;          // on phi
    };

    struct bem_element
    {
        double r = 0.0, chord = 0.0, theta = 0.0;   // radius, chord, twist + pitch [rad]
        double R = 1.0, R_hub = 0.0;                // tip and hub radius
        int B = 3;
        const airfoil *af = nullptr;
        double re = 0.0;                             // Reynolds number for the table blending (0: first table)
    };

    struct bem_result
    {
        double phi = 0.0, alpha = 0.0, a = 0.0, ap = 0.0, F = 1.0, W = 0.0;
        double cl = 0.0, cd = 0.0, cm = 0.0, cn = 0.0, ct = 0.0;
        double fn = 0.0, ft = 0.0, mz = 0.0;         // per unit span
        int region = 0;     // 0 no induction, 1 windmill (0, pi/2], 2 propeller brake, 3 (pi/2, pi)
        bool converged = true;
    };

    bem_result bem_solve(const bem_element&, double Vx, double Vy, double rho, const bem_options& = bem_options(), double skew=0.0);
    // the residual of Ning's method at phi (tests)
    double bem_residual(const bem_element&, double Vx, double Vy, double phi, const bem_options&, double *a=nullptr, double *ap=nullptr);

    struct blade_station
    {
        double r = 0.0;          // radius from the rotor axis
        double chord = 0.0;
        double twist = 0.0;      // rad
        int airfoil = 0;         // index into rotor_def::airfoils
    };

    struct rotor_def
    {
        int B = 3;
        double R_hub = 0.0, R_tip = 1.0;
        std::vector<blade_station> stations;     // increasing r
        std::vector<airfoil> airfoils;
    };

    struct rotor_result
    {
        double T = 0.0, Q = 0.0, P = 0.0;        // thrust, torque, power
        double Cp = 0.0, Ct = 0.0, lambda = 0.0;
        std::vector<bem_result> el;
    };

    // AeroDyn 15 blade file: NumBlNds rows BlSpn BlCrvAC BlSwpAC BlCrvAng BlTwist BlChord BlAFID [...]; the stations at
    // r = R_hub + BlSpn (the span from the blade root), twist in rad, airfoil = BlAFID - 1 (the order of AeroDyn's AFNames)
    bool read_aerodyn_blade(const std::string &file, double R_hub, std::vector<blade_station>&, std::string &msg);
    bool read_aerodyn_blade_string(const std::string &text, double R_hub, std::vector<blade_station>&, std::string &msg);

    rotor_result rotor_steady(const rotor_def&, double U, double Omega, double pitch, double rho, const bem_options& = bem_options());
}

#endif
