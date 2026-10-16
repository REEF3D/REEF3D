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

#include"coral/core/coral_gravity.h"
#include"coral/core/coral_model.h"

namespace coral
{

void gravity_load::eval(const model &m, double, VecX &f) const
{
    for(int i=0; i<m.nbody(); ++i)
        f.segment<3>(6*i) += m.bodies[i].mass*m.gravity;
}

double gravity_load::potential(const model &m) const
{
    double e = 0.0;
    for(const auto &b : m.bodies)
        e -= b.mass*m.gravity.dot(b.x);
    return e;
}

}
