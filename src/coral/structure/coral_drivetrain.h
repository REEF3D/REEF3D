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
// REEF3D::CORAL: drivetrain of a wind turbine: the rotor (hub) turns in the shaft joint (revolute, nacelle -> hub),
// the generator either on the same body (rigid drivetrain, its inertia N^2 J_gen added to the hub) or on a body of
// its own in a coaxial revolute joint (nacelle -> generator) with the shaft's torsion spring and damper between them:
//   T_shaft = k (theta_hub - theta_gen) + c (omega_hub - omega_gen)      (low-speed shaft)
// Generator torque on the high-speed shaft, T_gen(omega_hss), omega_hss = N omega (omega of the generator side
// relative to the nacelle); the low-speed shaft takes N T_gen against the rotation, the nacelle the reaction:
//   none          free
//   constant T    T_gen = T
//   slip          T_gen = slope (omega_hss - omega_sync)        (induction generator; a large slope holds the speed)
//   k_omega2      T_gen = K omega_hss^2 (omega_hss > 0), at most T_max   (region 2 of a variable-speed turbine)
// The torque can be set from outside each step (controllers, 3.5: set_torque_command). Electric power
// P = eta T_gen omega_hss.

#ifndef CORAL_DRIVETRAIN_H_
#define CORAL_DRIVETRAIN_H_

#include"coral/core/coral_load.h"
#include<string>

namespace coral
{
    class joint;

    enum class generator_law {none, constant, slip, k_omega2, command};

    struct drivetrain_settings
    {
        double N = 1.0;                  // gearbox ratio (high-speed / low-speed shaft)
        double k = 0.0, c = 0.0;         // torsion of the low-speed shaft (flexible drivetrain only)
        generator_law law = generator_law::none;
        double T = 0.0;                  // constant torque, or the cap of k_omega2 (0: none)
        double slope = 0.0, omega_sync = 0.0;   // slip: N m s / rad and rad/s on the high-speed shaft
        double K = 0.0;                  // k_omega2 [N m s^2 / rad^2]
        double eta = 1.0;                // generator efficiency
    };

    class drivetrain : public load
    {
    public:
        // gen: the generator's joint (nacelle -> generator body) of a flexible drivetrain, or nullptr (rigid)
        drivetrain(const joint &shaft, const joint *gen, const drivetrain_settings&, const std::string &name="drivetrain");

        std::string name() const override {return name_;}
        coupling_kind coupling() const override {return coupling_kind::local;}
        void coupled_bodies(std::vector<int>&) const override;
        void eval(const model&, double t, VecX &f) const override;
        void channel_names(std::vector<std::string>&) const override;
        void channel_values(const model&, double t, std::vector<double>&) const override;

        double omega_rotor(const model&) const;     // rad/s, hub relative to the nacelle
        double omega_hss(const model&) const;       // rad/s, generator side times N
        double generator_torque(const model&) const;
        double shaft_torque(const model&) const;
        void set_torque_command(double T) {cmd_ = T;}

        drivetrain_settings settings;

    private:
        std::string name_;
        const joint *shaft_, *gen_;
        double cmd_ = 0.0;
    };
}

#endif
