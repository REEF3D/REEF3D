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
// REEF3D::CORAL: wind turbine rotor in the multibody model (Phase 3.3): the hub turns in a revolute shaft joint on the
// nacelle (or the ground); the drivetrain (gearbox, generator inertia, optionally the shaft's torsion with the
// generator on a body of its own) and the generator torque law (coral_drivetrain.h); the aerodynamic load from blade
// element momentum at the stations of the blade definition (coral_rotor_aero.h). The blades:
//   rigid, on the hub: the blades' distributed mass as point masses of the hub along their span; collective pitch
//       of the aerodynamics (constant, a table, or set by a controller);
//   rigid or flexible, on pitch bearings (pitch_joints): one body per blade root in a revolute joint on the hub, turned
//       by a second-order actuator (coral_pitch_actuator.h, natural frequency and damping ratio), the pitch then in
//       the geometry; the bearing carries the blade's pitch inertia sum m (0.04^2 + 0.2^2) c^2 (c the chord; with
//       flexible blades the beams' torsional inertia adds as much, and the actuator is tuned to twice it);
//   flexible (blade_elements n): n co-rotational beam elements per blade from the root (the hub or the pitch bearing)
//       to the tip, with the mass and the flap and edge stiffness of the structural table (ElastoDyn's blade file)
//       oriented by its structural twist; EA = ea_per_mass * mass/length, GJ = gj_factor * EI_flap, rotary inertia per
//       length m (0.2 c)^2 about the chord normal and m (0.04 c)^2 about the chord (c the chord; ElastoDyn's blades
//       are rigid in extension and torsion), stiffness-proportional damping blade_beta; the aerodynamic stations on
//       the nearest node.
//
// Geometry: hub centre H, shaft axis e pointing downwind (tilt included by the user), the rotor turns positively
// about e (clockwise looking downwind) or the other way (ccw). Blade 1 points up (the vertical projected on the
// rotor plane) at azimuth 0, rotated by azimuth0; blade k follows at 2 pi k / B in the direction of rotation. Cone:
// the span direction cos(precone) r_k + sin(precone) e (positive: tips downwind; ElastoDyn's PreCone). The blade
// stations at the distance r from H along the span (r = R_hub + BlSpn of AeroDyn), twist and chord from the
// definition, pitch collective (constant, a table in time, or set by a controller).

#ifndef CORAL_ROTOR_H_
#define CORAL_ROTOR_H_

#include"coral/aero/coral_rotor_aero.h"
#include"coral/structure/coral_blade_table.h"
#include"coral/structure/coral_drivetrain.h"
#include"coral/structure/coral_pitch_actuator.h"
#include<memory>
#include<string>
#include<vector>

namespace coral
{
    class model;
    class joint;

    struct rotor_spec
    {
        std::string name = "rotor";
        int B = 3;
        Vec3 hub = Vec3::Zero();
        Vec3 axis = Vec3::UnitX();
        bool ccw = false;
        double precone = 0.0, azimuth0 = 0.0;          // rad
        double R_hub = 0.0;                            // distance of the blade roots from H along the span
        std::vector<blade_station> stations;           // r from H along the span
        std::vector<airfoil> airfoils;
        blade_structure structure;                     // empty: no blade mass
        double hub_mass = 0.0, hub_J = 0.0, hub_Jt = -1.0;   // hub inertia about the shaft and across it (-1: hub_J/2)
        int nacelle = -1;
        double omega0 = 0.0;                           // initial rotor speed [rad/s]
        double pitch = 0.0;                            // rad
        std::vector<double> pitch_t, pitch_v;          // collective pitch in time [rad]
        drivetrain_settings dt;
        double J_gen = 0.0;                            // generator inertia on the high-speed shaft
        bool pitch_joints = false;
        double pitch_freq = 2.0, pitch_zeta = 0.7;     // pitch actuator [Hz], damping ratio
        int blade_elements = 0;                        // > 0: flexible blades of this many beam elements
        double blade_beta = 0.0022;                    // stiffness-proportional damping of the blades [s]
        double ea_per_mass = 2.0e7, gj_factor = 0.3;
        rotor_aero_settings aero;                      // B, R_hub, R_tip set by the builder; bem, skew, tower, rho
    };

    struct rotor_handles
    {
        int hub = -1, gen = -1;
        joint *shaft = nullptr, *gen_joint = nullptr;
        rotor_aero *aero = nullptr;
        drivetrain *drive = nullptr;
        pitch_actuator *pitch = nullptr;
        std::vector<int> pitch_bodies;
        std::vector<joint*> pitch_joints;
        std::vector<std::vector<int>> blade_nodes;     // flexible blades: the node bodies from the root to the tip
        double R_tip = 0.0;                            // in the rotor plane
        double mass = 0.0, J = 0.0;                    // rotor mass (hub + blades) and inertia about the shaft (without the generator)
    };

    // adds the rotor to the model (wind: the inflow of the aerodynamic load)
    rotor_handles build_rotor(model&, std::shared_ptr<const wind_inflow>, const rotor_spec&);
}

#endif
