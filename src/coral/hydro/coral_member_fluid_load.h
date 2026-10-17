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
// REEF3D::CORAL: distributed loads of a resolved structure on a frame or superelement member (one-way,
// fluid_member in coral.txt). The fluid solver keeps the structure fixed (or prescribed) in its grid and
// integrates the pressure (and viscous) loads on its hull; the hull forces are binned in strips along the
// member axis by the projection of their points on it (ns strips of equal length between the member
// ends at the start, points beyond the ends to the end strips). Per strip the force goes to the member
// at the strip centre, split linearly between the two nodes of the element there.
// In time: samples at the fluid time levels, a Lagrange polynomial through the last order + 1 of them
// (as fluid_load); excitation: without the loads of the first sample (still water).

#ifndef CORAL_MEMBER_FLUID_LOAD_H_
#define CORAL_MEMBER_FLUID_LOAD_H_

#include"coral/core/coral_load.h"
#include"coral/core/coral_points.h"
#include<memory>
#include<string>
#include<vector>

namespace coral
{
    class member_fluid_load : public load
    {
    public:
        // elements: pairs of point indices of pts along the member from p1 to p2 (positions now)
        member_fluid_load(const model&, std::shared_ptr<const point_set> pts, const std::vector<std::pair<int,int>> &elements,
                          const Vec3 &p1, const Vec3 &p2, int nstrips, const std::string &name="fluid_member");

        std::string name() const override {return name_;}
        coupling_kind coupling() const override {return coupling_kind::local;}
        void coupled_bodies(std::vector<int> &b) const override {b.insert(b.end(), bodies_.begin(), bodies_.end());}
        void coupled_coordinates(std::vector<int> &c) const override {c.insert(c.end(), coords_.begin(), coords_.end());}
        void eval(const model&, double t, VecX &f) const override;
        void channel_names(std::vector<std::string>&) const override;
        void channel_values(const model&, double t, std::vector<double>&) const override;

        // the strips: axis from p1 (unit e, length L), their number
        const Vec3 &p1() const {return a_;}
        const Vec3 &axis() const {return e_;}
        double length() const {return L_;}
        int nstrips() const {return ns_;}

        // strip forces (3 per strip, global) at the fluid time t
        void set_sample(double t, const VecX &F);
        bool has_sample() const {return !ts_.empty();}
        VecX forces(double t) const;
        // total force and moment about the origin (global) at time t
        void total(const model&, double t, Vec3 &F, Vec3 &M) const;

        int order = 2;
        bool excitation = false;

    private:
        std::string name_;
        std::shared_ptr<const point_set> pts_;
        Vec3 a_, e_;
        double L_;
        int ns_;
        std::vector<int> pa_, pb_;      // per strip: the two points of its element
        std::vector<double> w_;         // weight of pb
        std::vector<int> bodies_, coords_;
        std::vector<double> ts_;
        std::vector<VecX> Fs_;
        VecX F0_;
    };
}

#endif
