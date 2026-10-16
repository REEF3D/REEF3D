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
// REEF3D::CORAL: Morison loads distributed along flexible members (frame members on node bodies,
// superelement members on recovered node motion), strip theory per element between two material points
// of a point_set. Per unit length, normal to the current element axis:
//   f = rho A Cm a_w - rho A Ca a_s + 1/2 rho Cd D |u_r| u_r,   A = pi D^2/4,  u_r = u_w + U_c - u_s
// plus the buoyancy rho g A (upwards) of the submerged strips. Positions, velocities and accelerations of
// a strip are interpolated linearly between the element ends; its force goes back to the two ends with
// the same weights. Wave kinematics from the wave_field (linear, z = 0 kinematics above the still water
// level); only the parts of the strips below the instantaneous elevation carry load (a strip cut by the
// surface counts with its wet part, so the load is continuous in time). Axial (end) effects are left out.

#ifndef CORAL_MEMBER_MORISON_H_
#define CORAL_MEMBER_MORISON_H_

#include"coral/core/coral_load.h"
#include"coral/core/coral_points.h"
#include"coral/hydro/coral_waves.h"
#include<memory>
#include<string>
#include<utility>
#include<vector>

namespace coral
{
    struct member_morison_coef
    {
        double D = 1.0;
        double Cd = 1.0;
        double Ca = 1.0;
        double Cm = -1.0;          // < 0: 1 + Ca
        int n = 2;                 // strips per element
        bool buoyancy = true;
    };

    class member_morison : public load
    {
    public:
        // elements: pairs of point indices of pts
        member_morison(std::shared_ptr<const point_set> pts, const std::vector<std::pair<int,int>> &elements,
                       const member_morison_coef&, std::shared_ptr<const wave_field>, double rho, double g,
                       const Vec3 &current, const std::string &name);

        std::string name() const override {return name_;}
        coupling_kind coupling() const override {return coupling_kind::local;}
        void coupled_bodies(std::vector<int> &b) const override {b.insert(b.end(), bodies_.begin(), bodies_.end());}
        void coupled_coordinates(std::vector<int> &c) const override {c.insert(c.end(), coords_.begin(), coords_.end());}
        void eval(const model&, double t, VecX &f) const override;
        void channel_names(std::vector<std::string>&) const override;
        void channel_values(const model&, double t, std::vector<double>&) const override;

        // total force (global) and moment about the origin of all strips
        void total(const model&, double t, Vec3 &F, Vec3 &M) const;

    private:
        template<class Sink> void strips(const model&, double t, Sink sink) const;
        std::string name_;
        std::shared_ptr<const point_set> pts_;
        std::vector<std::pair<int,int>> elem_;
        std::vector<int> pidx_;        // the two end points of each element
        member_morison_coef c_;
        std::shared_ptr<const wave_field> waves_;
        double rho_, g_;
        Vec3 current_;
        std::vector<int> bodies_, coords_;
    };
}

#endif
