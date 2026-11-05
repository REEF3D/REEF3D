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
// REEF3D::CORAL: blade pitch actuators of a wind turbine: one revolute pitch joint per blade (hub -> pitch bearing,
// the axis such that a positive rotation turns the blade towards feather), each driven to its command by a second-order
// actuator, a torque about the joint
//   T = -k (theta - theta_cmd) - c (theta' - theta_cmd'),   theta = theta0 + q_joint
// with k = J omega^2, c = 2 zeta J omega (J: the blade's inertia about the pitch axis, omega = 2 pi f the actuator's
// natural frequency). The command: a constant, a collective table in time, or set per blade (controllers, 3.5).

#ifndef CORAL_PITCH_ACTUATOR_H_
#define CORAL_PITCH_ACTUATOR_H_

#include"coral/core/coral_load.h"
#include<string>
#include<vector>

namespace coral
{
    class joint;

    class pitch_actuator : public load
    {
    public:
        pitch_actuator(std::vector<const joint*> joints, double k, double c, double theta0, const std::string &name="pitch");

        std::string name() const override {return name_;}
        coupling_kind coupling() const override {return coupling_kind::local;}
        void coupled_bodies(std::vector<int>&) const override;
        void eval(const model&, double t, VecX &f) const override;
        void channel_names(std::vector<std::string>&) const override;
        void channel_values(const model&, double t, std::vector<double>&) const override;

        void set_command(int blade, double theta) {cmd_[blade] = theta; table_ = false;}
        void set_command_all(double theta) {for(double &c : cmd_) c = theta; table_ = false;}
        void set_table(std::vector<double> t, std::vector<double> th) {tt_ = std::move(t); tv_ = std::move(th); table_ = !tt_.empty();}
        // command and its rate at t
        void command(int blade, double t, double &th, double &rate) const;
        double angle(const model&, int blade) const;
        double rate(const model&, int blade) const;
        double torque(const model&, int blade, double t) const;
        int blades() const {return (int)j_.size();}

        double k, c;

    private:
        std::vector<const joint*> j_;
        double theta0_;
        std::vector<double> cmd_, tt_, tv_;
        bool table_ = false;
        std::string name_;
    };
}

#endif
