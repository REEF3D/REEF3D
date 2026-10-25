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
// REEF3D::CORAL: distributed structural properties of a blade along its span (from the root to the tip): fraction of
// the blade length, pitch axis (fraction of the chord from the leading edge), structural twist, mass per length,
// flapwise and edgewise bending stiffness. Readers: ElastoDyn's blade file (NBlInpSt rows BlFract PitchAxis StrcTwst
// BMassDen FlpStff EdgStff, the factors AdjBlMs, AdjFlSt, AdjEdSt applied) or a plain table of these six columns.

#ifndef CORAL_BLADE_TABLE_H_
#define CORAL_BLADE_TABLE_H_

#include<string>
#include<vector>

namespace coral
{
    struct blade_structure
    {
        std::vector<double> frac, pitch_axis, twist, mass, EI_flap, EI_edge;   // twist in rad

        // linear in the fraction (held at the ends)
        double mass_at(double f) const;
        // mass of the blade of length L
        double total_mass(double L) const;

        static bool read(const std::string &file, blade_structure&, std::string &msg);
        static bool read_string(const std::string &text, blade_structure&, std::string &msg);
    };
}

#endif
