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
// REEF3D::CORAL: linear 6x6 stiffness and damping of a body against the ground, about a reference
// point P (e.g. a linearized mooring system or a foundation matrix):
//   F = F0 - K xi - C dxi/dt,  xi = [p - p0; theta],  theta = log(R R0^T) (global rotation vector)
// F = [force (global) at the current reference point; moment (global)].

#ifndef CORAL_LINEAR6_H_
#define CORAL_LINEAR6_H_

#include"coral/core/coral_load.h"
#include"coral/core/coral_types.h"

namespace coral
{
    class linear6 : public load
    {
    public:
        linear6(const model&, int body, const Vec3 &P, const Mat6 &K, const Mat6 &C,
                const Vec6 &F0=Vec6::Zero(), const std::string &name="linear6");

        std::string name() const override {return name_;}
        void eval(const model&, double t, VecX &f) const override;
        double potential(const model&) const override;

        Vec6 displacement(const model&) const;
        Vec6 velocity(const model&) const;

    private:
        std::string name_;
        int b_;
        Vec3 s_, p0_;
        Mat3 R0_;
        Mat6 K_, C_;
        Vec6 F0_;
    };
}

#endif
