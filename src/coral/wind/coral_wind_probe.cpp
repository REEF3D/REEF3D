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

#include"coral/wind/coral_wind_probe.h"
#include"coral/core/coral_model.h"

namespace coral
{

wind_probe::wind_probe(const model &m, std::shared_ptr<const wind_inflow> w, int b, const Vec3 &p, const std::string &name)
    : wind_(std::move(w)), b_(b), s_(m.body(b).to_body(p)), name_(name) {}

void wind_probe::channel_names(std::vector<std::string> &c) const
{
    for(const char *x : {"u", "v", "w", "U"})
        c.push_back(name_ + "_" + x);
}

void wind_probe::channel_values(const model &m, double t, std::vector<double> &v) const
{
    const Vec3 u = wind_->velocity(m.body(b_).point(s_), t);
    v.push_back(u(0)); v.push_back(u(1)); v.push_back(u(2)); v.push_back(u.norm());
}

}
