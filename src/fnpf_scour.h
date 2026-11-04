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
Author: Hans Bihs
--------------------------------------------------------------------*/

#ifndef FNPF_SCOUR_H_
#define FNPF_SCOUR_H_

#include "scour_core.h"

class fdm_fnpf;

// P 187, P 198: the scour parameters of FNPF (scour_core.h). The bed (slip) velocity of the potential
// flow is the undisturbed near-bed orbital velocity the engineering formulas expect; there is no
// computed current (P 189 or P 196 give the ambient current).
class fnpf_scour : public scour_core
{
public:
    fnpf_scour(lexer *, fdm_fnpf *, ghostcell *);
    virtual ~fnpf_scour() = default;

    // every time step: the bed velocity and its statistics
    void sample(lexer *, fdm_fnpf *, ghostcell *);

private:
    void bedvel(lexer *, fdm_fnpf *);
};

#endif
