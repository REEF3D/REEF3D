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
// REEF3D::CORAL: prescribed motion of a joint coordinate (forced oscillation of a body in a resolved fluid:
// added mass and damping, the hydrodynamic database from REEF3D runs). A stiff spring-damper in the joint
// drives the coordinate to
//   q_p(t) = q0 + a r(t) sin(2 pi t / T + phi),   r(t) = (1 - cos(pi t / t_r))/2 for t < t_r, else 1
//   F = -k (q - q_p) - c (q' - q_p')
// (force along the axis of a prismatic or cylindrical joint, torque about the axis of a revolute one).
// With k much larger than the stiffness of the rest of the system the motion error is F/k; F is then the
// force needed to impose the motion (the joint reaction along the axis).

#ifndef CORAL_JOINT_DRIVE_H_
#define CORAL_JOINT_DRIVE_H_

#include"coral/core/coral_load.h"
#include<string>

namespace coral
{
    class joint;

    struct drive_settings
    {
        double q0 = 0.0, amp = 0.0, period = 1.0, phase = 0.0, ramp = 0.0;
        double k = 1.0e7, c = 1.0e4;
    };

    class joint_drive : public load
    {
    public:
        // the joint must live as long as the load (it is owned by the model)
        joint_drive(const joint&, bool rotational, const drive_settings&, const std::string &name="drive");

        std::string name() const override {return name_;}
        coupling_kind coupling() const override {return coupling_kind::local;}
        void coupled_bodies(std::vector<int>&) const override;
        void eval(const model&, double t, VecX &f) const override;
        void channel_names(std::vector<std::string>&) const override;
        void channel_values(const model&, double t, std::vector<double>&) const override;

        // prescribed coordinate, rate and acceleration at t
        void prescribed(double t, double &q, double &qd, double &qdd) const;
        double force(const model&, double t) const;

        drive_settings settings;

    private:
        std::string name_;
        const joint *j_;
        bool rot_;
    };
}

#endif
