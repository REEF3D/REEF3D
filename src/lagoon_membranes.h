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

#ifndef LAGOON_MEMBRANES_H_
#define LAGOON_MEMBRANES_H_

// The membranes of a run (X 330, NHFLOW) in its LAGOON store (P 18; LAGOON SPEC.md
// section 10): one particle set of polygons, particles/nhflow_membrane with role
// "membrane", holding the membranes of every output side by side as in the
// REEF3D-NHFLOW-Membrane .vtp files - the nodes as its points with velocity,
// displacement and load (per area), the triangles as its cells with the number of
// their membrane ("membrane"), dp, force, tension and tag.
//
// Every membrane prints for itself (rank 0). An output goes into the store once a
// membrane prints in a later time step, or at the end of the run (finish). A
// membrane that did not print in a stored step (its own 'print dt') is stored with
// its last printed state, so no membrane leaves the picture.

#include<cstdint>
#include<vector>

class lexer;

class lagoon_membranes
{
public:
    struct panels
    {
        std::vector<float> xyz;           // 3 per node
        std::vector<float> velocity;      // 3 per node
        std::vector<float> displacement;  // 3 per node
        std::vector<float> load;          // 3 per node, per area
        std::vector<int32_t> triangles;   // 3 node numbers per triangle
        std::vector<float> dp;            // 1 per triangle
        std::vector<float> force;         // 3 per triangle
        std::vector<float> tension;       // 1 per triangle
        std::vector<int32_t> tag;         // 1 per triangle
    };

    // rank 0: membrane `membrane` printed (its .vtp file, number num). True when the
    // output goes into the store (with P 18 1 the .vtp file is then left out);
    // false without P 18 or once the store could not be written.
    static bool output(lexer*, int membrane, int num, panels &&m);

    // rank 0, at the end of the run (lagoon_output::finish_all): the last output
    static void finish(lexer*);
};

#endif
