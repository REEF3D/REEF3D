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
// REEF3D::CORAL: constant point force F (global) at a point of a body and a constant moment M (global),
// optionally switched on linearly over a ramp time (static pushover, prescribed loads, tests).

#ifndef CORAL_FORCE_H_
#define CORAL_FORCE_H_

#include"coral/core/coral_load.h"
#include"coral/core/coral_types.h"
#include<string>

namespace coral
{
    class model;

    class point_force : public load
    {
    public:
        // p: global point on body b at creation
        point_force(const model&, int b, const Vec3 &p, const Vec3 &F, const Vec3 &M, double ramp=0.0, const std::string &name="force");

        std::string name() const override {return name_;}
        coupling_kind coupling() const override {return coupling_kind::local;}
        void coupled_bodies(std::vector<int> &b) const override {b.push_back(b_);}
        void eval(const model&, double t, VecX &f) const override;

        Vec3 F, M;

    private:
        std::string name_;
        int b_;
        Vec3 s_;
        double ramp_;
    };
}

#endif
