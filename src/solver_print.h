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

#ifndef SOLVER_PRINT_H_
#define SOLVER_PRINT_H_

// The line of the solvers per solve (REEFMG cycles, Fi_iter) off while the 6DOF bodies of FNPF solve their psi
// problems (fnpf_6DOF::solve_psi): with many bodies they were most of the terminal output; their number and
// iterations come as one line per step instead (fnpf_6DOF::psi_report)
inline bool &solver_print_quiet()
{
    static bool quiet = false;
    return quiet;
}

#endif
