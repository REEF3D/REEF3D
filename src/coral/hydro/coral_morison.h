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
// REEF3D::CORAL: Morison loads on slender cylindrical members of a body (strip theory). Per unit length,
// normal to the member axis:
//   f = rho A Cm a_w - rho A Ca a_b + 1/2 rho Cd D |u_r| u_r,   A = pi D^2/4,  u_r = u_w + U_c - u_b
// a_w, u_w wave acceleration and velocity (wave_field, linear, z = 0 kinematics above the still
// water level), U_c a uniform current, a_b, u_b of the member. Only the part below the instantaneous
// wave elevation carries load. The member is integrated with n strips (midpoint rule). Members on the
// ground (body -1) give forces through total() only (fixed structures without dynamics).

#ifndef CORAL_MORISON_H_
#define CORAL_MORISON_H_

#include"coral/core/coral_load.h"
#include"coral/hydro/coral_waves.h"
#include<memory>
#include<string>
#include<vector>

namespace coral
{
    struct morison_member
    {
        Vec3 p1 = Vec3::Zero(), p2 = Vec3::Zero();    // end points (global, at creation)
        double D = 1.0;
        double Cd = 1.0;
        double Ca = 1.0;
        double Cm = 2.0;      // 1 + Ca for a member in waves (Froude-Krylov + diffraction)
        int n = 20;
    };

    class morison_load : public load
    {
    public:
        morison_load(const model&, int body, const std::vector<morison_member>&, std::shared_ptr<const wave_field>,
                     double rho=1025.0, const Vec3 &current=Vec3::Zero(), const std::string &name="morison");

        std::string name() const override {return name_;}
        coupling_kind coupling() const override {return coupling_kind::local;}
        void coupled_bodies(std::vector<int> &b) const override {b.push_back(b_);}
        void eval(const model&, double t, VecX &f) const override;
        void channel_names(std::vector<std::string>&) const override;
        void channel_values(const model&, double t, std::vector<double>&) const override;

        // total force (global) and moment about the origin (global) of the members
        void total(const model&, double t, Vec3 &F, Vec3 &M) const;

    private:
        struct strip {Vec3 s; Vec3 e; double L; double D, Cd, Ca, Cm;};    // body frame: centre, axis
        bool strip_force(const model&, const strip&, double t, Vec3 &p, Vec3 &f) const;
        std::string name_;
        int b_;
        std::vector<strip> strips_;
        std::shared_ptr<const wave_field> waves_;
        double rho_;
        Vec3 current_;
    };
}

#endif
