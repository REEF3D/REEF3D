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
// REEF3D::CORAL: wind probe: the velocity of the inflow at a point fixed in space or on a body, as output
// channels <name>_u, _v, _w (global axes) and <name>_U (magnitude). It adds no force.

#ifndef CORAL_WIND_PROBE_H_
#define CORAL_WIND_PROBE_H_

#include"coral/core/coral_load.h"
#include"coral/wind/coral_wind.h"
#include<memory>
#include<string>

namespace coral
{
    class wind_probe : public load
    {
    public:
        // p: global point (on body b at the current configuration, or fixed for b = -1)
        wind_probe(const model&, std::shared_ptr<const wind_inflow>, int b, const Vec3 &p, const std::string &name);

        std::string name() const override {return name_;}
        coupling_kind coupling() const override {return coupling_kind::local;}
        void eval(const model&, double, VecX&) const override {}
        void channel_names(std::vector<std::string>&) const override;
        void channel_values(const model&, double t, std::vector<double>&) const override;

    private:
        std::shared_ptr<const wind_inflow> wind_;
        int b_;
        Vec3 s_;
        std::string name_;
    };
}

#endif
