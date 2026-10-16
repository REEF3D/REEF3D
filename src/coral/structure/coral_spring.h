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
// REEF3D::CORAL: spring-damper between two points (body a, body b; -1 ground), e.g. a tether.
//   tension T = T0 + k (L - L0) + c dL/dt,  L = |p_b - p_a|
//   tension_only: T = max(T, 0) (slack line)
// Points are given in the global frame in the configuration when the spring is created.
// L0 < 0: the length at creation.

#ifndef CORAL_SPRING_H_
#define CORAL_SPRING_H_

#include"coral/core/coral_load.h"
#include"coral/core/coral_types.h"

namespace coral
{
    class spring_damper : public load
    {
    public:
        spring_damper(const model&, int a, const Vec3 &pa, int b, const Vec3 &pb,
                      double k, double c, double L0=-1.0, double T0=0.0, bool tension_only=false,
                      const std::string &name="spring");

        std::string name() const override {return name_;}
        coupling_kind coupling() const override {return coupling_kind::local;}
        void coupled_bodies(std::vector<int> &b) const override {b.push_back(a_); b.push_back(b_);}
        void eval(const model&, double t, VecX &f) const override;
        double potential(const model&) const override;
        void channel_names(std::vector<std::string>&) const override;
        void channel_values(const model&, double t, std::vector<double>&) const override;

        double length(const model&) const;
        double tension(const model&) const;

    private:
        std::string name_;
        int a_, b_;
        Vec3 sa_, sb_;
        double k_, c_, L0_, T0_;
        bool tension_only_;
    };
}

#endif
