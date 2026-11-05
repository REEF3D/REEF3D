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
// REEF3D::CORAL: output of the linear analysis. modes.dat: one line per mode with the frequency, the
// period, the damping ratio and the split of the modal kinetic energy over the six global directions
// (translation x y z, rotation about x y z), so that fore-aft, side-side, torsion and heave modes can
// be told apart. mode_shapes.dat: per mode one line per body with its name, its position and the six
// components of the shape (translation global, rotation global).

#ifndef CORAL_MODES_OUT_H_
#define CORAL_MODES_OUT_H_

#include"coral/core/coral_types.h"
#include<array>
#include<string>

namespace coral
{
    class model;
    struct mode_result;

    // fractions of the modal kinetic energy of mode q: x y z rx ry rz (sum 1)
    std::array<double,6> energy_split(const model&, const mode_result&, int q);

    // short name of the dominant direction ("x", "y", "z", "rx", "ry", "rz")
    std::string dominant(const std::array<double,6>&);

    // modes.dat and mode_shapes.dat in the folder dir (with the trailing '/' or empty)
    void write_modes(const std::string &dir, const model&, const mode_result&);

    // the table on stdout
    void print_modes(const model&, const mode_result&, bool damped);
}

#endif
