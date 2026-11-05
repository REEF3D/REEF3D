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
// REEF3D::CORAL: prescribed force (or torque) along the axis of a joint, from a table t F (system
// identification of a WEC: a multisine or chirp force of the PTO motor, e.g. the intrinsic admittance
// Y_i = V / F of the relative motion in calm water):
//   F(t) = s r(t) F_tab(t),   r(t) = (1 - cos(pi t / t_r))/2 for t < t_r, else 1
// F_tab is interpolated linearly in the table and 0 outside it (periodic: the table repeats with the
// period t_last - t_first). The force acts like a PTO force: along the axis of a prismatic or
// cylindrical joint on body b and against it on body a, a torque about the axis of a revolute one.

#ifndef CORAL_JOINT_FORCE_H_
#define CORAL_JOINT_FORCE_H_

#include"coral/core/coral_load.h"
#include<string>
#include<vector>

namespace coral
{
    class joint;

    struct jforce_settings
    {
        std::vector<double> t, F;       // table, t increasing
        double scale = 1.0, ramp = 0.0;
        bool periodic = false;
    };

    class joint_force : public load
    {
    public:
        // the joint must live as long as the load (it is owned by the model)
        joint_force(const joint&, bool rotational, const jforce_settings&, const std::string &name="jforce");

        std::string name() const override {return name_;}
        coupling_kind coupling() const override {return coupling_kind::local;}
        void coupled_bodies(std::vector<int>&) const override;
        void eval(const model&, double t, VecX &f) const override;
        void channel_names(std::vector<std::string>&) const override;
        void channel_values(const model&, double t, std::vector<double>&) const override;

        // the prescribed force at t
        double force(double t) const;
        // reads a table "t F" (blank lines and lines from '#' on are skipped); false with msg on an error
        static bool read_table(const std::string &file, std::vector<double> &t, std::vector<double> &F, std::string &msg);

        jforce_settings settings;

    private:
        std::string name_;
        const joint *j_;
        bool rot_;
    };
}

#endif
