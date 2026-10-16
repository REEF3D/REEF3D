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
// REEF3D::CORAL: builder of beam structures (frames): members between two points, divided into n
// co-rotational beam elements (coral_beam) on node bodies. A node is a small rigid body (6 DOF) named
// "<member>:<i>" (i = 0 .. n); members that meet at a point share the node. The mass of the elements
// is lumped on the nodes (half an element each, with the rotary inertia of the section and of the
// half element). A member end can also sit on a rigid body or the ground (no node there).

#ifndef CORAL_FRAME_H_
#define CORAL_FRAME_H_

#include"coral/core/coral_types.h"
#include"coral/structure/coral_beam.h"
#include<string>
#include<vector>

namespace coral
{
    class model;

    struct member_spec
    {
        std::string name;
        Vec3 p1 = Vec3::Zero(), p2 = Vec3::Zero();
        int n = 10;
        beam_section s1, s2;           // sections at p1 and p2 (linear in between)
        bool tapered = false;          // false: s1 everywhere
        Vec3 orient = Vec3::UnitZ();   // local y axis (projected normal to the member)
        double beta = 0.0;             // stiffness-proportional damping [s]
        int body1 = -2, body2 = -2;    // -2: node; -1: ground; >= 0: the end sits on that rigid body
    };

    class frame_builder
    {
    public:
        explicit frame_builder(model &m, double tol=1.0e-9) : m_(m), tol_(tol) {}

        // adds the member; returns the model body index of each station (-1/-2 for ends on the ground
        // or a body). The elements are added to the model as loads.
        std::vector<int> add_member(const member_spec&);

        // the node at p (existing within tol, else a new massless node named name)
        int node_at(const Vec3 &p, const std::string &name);

        // moves the lumped masses into the node bodies; call once after the last member
        void finish();

        const std::vector<int> &nodes() const {return nodes_;}

    private:
        model &m_;
        double tol_;
        std::vector<int> nodes_;
        std::vector<double> mass_;
        std::vector<Mat3> inertia_;
        bool finished_ = false;
    };
}

#endif
