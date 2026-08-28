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
Author: Elyas Larkermani
--------------------------------------------------------------------*/

// In-process Chrono contact for CFD 6DOF bodies (X 13 1).
// Chrono owns gravity-free rigid-body EOM and contact (NSC hard or SMC soft).
// REEF3D hydro wrenches (already include gravity) are applied each step.

#ifndef SIXDOF_CHRONO_H_
#define SIXDOF_CHRONO_H_

#include<vector>

class lexer;
class fdm;
class ghostcell;
class sixdof_obj;

class sixdof_chrono
{
public:
    sixdof_chrono();
    ~sixdof_chrono();

    void initialize(lexer*, std::vector<sixdof_obj*>&);
    void advance(lexer*, fdm*, ghostcell*, std::vector<sixdof_obj*>&, int, bool);

private:
    struct impl;
    impl *pimpl;
};

#endif
