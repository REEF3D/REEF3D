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
// REEF3D::CORAL: quasi-static elastic catenary of a mooring line with a seabed portion and seabed
// friction (Jonkman 2007, the formulation of MAP/FAST): anchor at the origin of the line's vertical
// plane, fairlead at horizontal distance xF and height zF above the anchor.
//   L unstretched length, w submerged weight per length, EA axial stiffness, CB seabed friction
// HF, VF horizontal and vertical force at the fairlead (tension components pulling the fairlead
// towards the anchor horizontally and down). A line with length to spare (xF + zF <= L with a seabed)
// is slack: HF = 0, it hangs vertically from the fairlead and the rest lies loosely on the seabed.

#ifndef CORAL_CATENARY_H_
#define CORAL_CATENARY_H_

namespace coral
{
    struct catenary_line
    {
        double L = 0.0, w = 0.0, EA = 0.0, CB = 0.0;
    };

    // fairlead position (xF, zF) for the fairlead forces (HF > 0, VF)
    void catenary_ends(const catenary_line&, double HF, double VF, double &xF, double &zF);

    // fairlead forces for the fairlead position; HF, VF are the start values on entry (<= 0: own
    // start values); false if Newton does not converge
    bool catenary_solve(const catenary_line&, double xF, double zF, double &HF, double &VF);

    // length on the seabed for the fairlead forces (0: line hanging free)
    double catenary_seabed_length(const catenary_line&, double HF, double VF);

    // point (x, z) at the unstretched arc length s from the anchor
    void catenary_profile(const catenary_line&, double HF, double VF, double s, double &x, double &z);

    // anchor forces
    void catenary_anchor(const catenary_line&, double HF, double VF, double &HA, double &VA);
}

#endif
