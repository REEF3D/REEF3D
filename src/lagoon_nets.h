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

#ifndef LAGOON_NETS_H_
#define LAGOON_NETS_H_

// The nets of a run in its LAGOON store (P 18; LAGOON SPEC.md section 10): the bar
// nets (X 320 types 1 - 3 and 12 - 14) of every output as one set of lines,
// particles/<solver>_net with role "net" - the knots as its points with Tension
// (and, for dynamic nets, Velocity and Acceleration) as in the REEF3D-Net .vtk
// files, the bars as its cells with the number of their net.
//
// Every net prints for itself (rank 0). An output goes into the store once the
// first net prints the next one, or at the end of the run (finish), so a net printed
// twice with the same number (the start: the initialisation, then the first step)
// is stored as its .vtk file is left, the later state.

#include<cstdint>
#include<vector>

class lexer;

class lagoon_nets
{
public:
    struct knots
    {
        std::vector<float> xyz;           // 3 per knot
        std::vector<int32_t> bars;        // 2 knot numbers per bar
        std::vector<float> tension;       // 1 per knot
        std::vector<float> velocity;      // 3 per knot, or empty (a quasi-static net)
        std::vector<float> acceleration;  // 3 per knot, or empty
    };

    // rank 0: net `net` printed output num (the number of its .vtk file). True when
    // the output goes into the store (with P 18 1 the .vtk file is then left out);
    // false without P 18 or once the store could not be written.
    static bool output(lexer*, int net, int num, knots &&k);

    // rank 0, at the end of the run (lagoon_output::finish_all): the last output
    static void finish(lexer*);
};

#endif
