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
// REEF3D::CORAL: quadratic (viscous) damping of a body, e.g. the drag a potential-flow solver (FNPF, a
// hydrodynamic database) leaves out:
//   F_i = -Q_i |u_i| u_i,   u = [velocity of the point P (global); angular velocity (global)]
// Q from the keyword directly (diagonal, surge sway heave roll pitch yaw) or Q = rho Cd A / 2 along an
// axis e (u_e = u.e, F = -rho Cd A |u_e| u_e e / 2, as for the drag on the waterplane area of a heaving
// float). The force acts at P (default: the CoG), the moment about the body.

#ifndef CORAL_BODY_DRAG_H_
#define CORAL_BODY_DRAG_H_

#include"coral/core/coral_load.h"
#include"coral/core/coral_types.h"
#include<string>

namespace coral
{
    class body_drag : public load
    {
    public:
        // Q: diagonal quadratic coefficients (N s^2/m^2 and N m s^2); P: point of the translational part
        body_drag(const model&, int body, const Vec3 &P, const Vec6 &Q, const std::string &name="drag");
        // along the axis e (unit vector, global): Q_e = rho Cd A / 2
        body_drag(const model&, int body, const Vec3 &P, const Vec3 &e, double Qe, const std::string &name="drag");

        std::string name() const override {return name_;}
        coupling_kind coupling() const override {return coupling_kind::local;}
        void coupled_bodies(std::vector<int> &b) const override {b.push_back(b_);}
        void eval(const model&, double t, VecX &f) const override;

        // the force (global) and moment (global) of the current state, and its power (negative)
        Vec6 force(const model&) const;
        double power(const model&) const;

        void channel_names(std::vector<std::string> &c) const override;
        void channel_values(const model&, double t, std::vector<double> &v) const override;

    private:
        std::string name_;
        int b_;
        Vec3 s_;
        Vec6 Q_ = Vec6::Zero();
        Vec3 e_ = Vec3::Zero();
        double Qe_ = 0.0;
        bool axis_ = false;
    };
}

#endif
