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
// REEF3D::CORAL: electrical side of a power take-off (wave-to-wire): the mechanical power of the PTO
// through a gear or motion conversion to the generator shaft and on to the electrical power at the
// grid side of the power electronics. The generator force is set by the PTO (its control law), the
// machine is current-controlled (the electrical dynamics are fast against the waves), so the model
// gives the power flow and its losses without feeding back on the motion.
//
//   gear / conversion   shaft speed w = ratio qd, generator torque T = -F/ratio g    (rack and pinion,
//                       ball screw, gearbox: ratio = 1/r_pinion, 2 pi/lead, n; a linear generator: 1)
//                       g = eta_gear generating, 1/eta_gear motoring, so P_shaft = g P_mech
//   machine             ideal: P_el = eta P_shaft (generating), P_shaft/eta (motoring)
//                       map:   the same with eta(|w|, |T|) from a table (bilinear, end values outside)
//                       pmsg:  permanent-magnet synchronous machine (rotary or linear) with field-oriented
//                              control (i_d = 0): i_q = T/kt, kt = 3/2 p lambda (rotary, p pole pairs) or
//                              3/2 pi/tau lambda (linear, pole pitch tau), P_el = P_shaft - 3/2 R i_q^2
//                              - k_fe w^2 (copper and iron losses; holding a force costs power)
//   converter           P_grid = eta_conv P_el (P_el > 0), P_el/eta_conv (P_el < 0)
// P_mech = -F qd is the absorbed power of the PTO, F the generator force on body b along +q.
// With reactive control the power flows back in parts of the cycle and pays the losses twice.

#ifndef CORAL_GENERATOR_H_
#define CORAL_GENERATOR_H_

#include<string>
#include<vector>

namespace coral
{
    enum class generator_kind {ideal, map, pmsg};

    struct generator_settings
    {
        generator_kind kind = generator_kind::ideal;
        double ratio = 1.0;          // shaft speed per unit PTO rate
        double eta_gear = 1.0;       // gear or conversion efficiency
        double eta = 1.0;            // ideal machine
        double R = 0.0;              // pmsg: phase resistance [Ohm]
        double kt = 0.0;             // pmsg: torque [N m/A] or force constant [N/A], 3/2 p lambda or 3/2 pi/tau lambda
        double k_fe = 0.0;           // pmsg: iron and eddy current losses k_fe w^2
        double I_max = 0.0;          // pmsg: current limit (for the channel only), 0: none
        double eta_conv = 1.0;       // power electronics
        // map: speeds (rows), torques (columns), efficiencies row by row
        std::vector<double> map_w, map_T, map_eta;
    };

    struct generator_power
    {
        double P_mech = 0.0;     // absorbed by the PTO
        double P_shaft = 0.0;    // at the generator shaft (after the gear)
        double P_el = 0.0;       // at the machine terminals
        double P_grid = 0.0;     // after the power electronics (the delivered power)
        double w = 0.0;          // shaft speed
        double T = 0.0;          // shaft torque of the generator (positive: brakes positive w)
        double I = 0.0;          // pmsg: q-axis current amplitude
    };

    class generator_model
    {
    public:
        generator_model() = default;
        explicit generator_model(const generator_settings &s) : settings(s) {}

        // F: generator force on the moving body along +q, qd: the PTO rate
        generator_power convert(double F, double qd) const;
        // efficiency of the map at |w|, |T|
        double map_eta(double w, double T) const;
        // reads the map: first line 0 T_1 .. T_m, then lines w_i eta_i1 .. eta_im ('#' comments);
        // false with msg on an error
        static bool read_map(const std::string &file, generator_settings &s, std::string &msg);
        // checks the settings, empty if fine
        std::string check() const;

        generator_settings settings;
    };
}

#endif
