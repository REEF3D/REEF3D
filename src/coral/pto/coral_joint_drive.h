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
// Instead of the sine, any motion from a table t q (file): q_p(t) = q0 + s r(t) Q(t), Q the C1 cubic
// Hermite curve through the table points (Catmull-Rom slopes), held at the end values outside the
// table unless periodic (repeats with the period t_last - t_first); measured motions, multisines,
// irregular forced motions (WEC-Sim's PTO actuation motion).
// With k much larger than the stiffness of the rest of the system the motion error is F/k; F is then the
// force needed to impose the motion (the joint reaction along the axis).

#ifndef CORAL_JOINT_DRIVE_H_
#define CORAL_JOINT_DRIVE_H_

#include"coral/core/coral_load.h"
#include<string>
#include<vector>

namespace coral
{
    class joint;

    struct drive_settings
    {
        double q0 = 0.0, amp = 0.0, period = 1.0, phase = 0.0, ramp = 0.0;
        double k = 1.0e7, c = 1.0e4;
        std::vector<double> t, q;        // table motion (empty: the sine)
        double scale = 1.0;
        bool periodic = false;
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
        // the table curve and its first and second derivatives at t (without q0, scale and ramp)
        void table(double t, double &Q, double &Qd, double &Qdd) const;
        // reads a table "t q" (blank lines and lines from '#' on are skipped); false with msg on an error
        static bool read_table(const std::string &file, std::vector<double> &t, std::vector<double> &q, std::string &msg);

        drive_settings settings;

    private:
        std::string name_;
        const joint *j_;
        bool rot_;
    };
}

#endif
