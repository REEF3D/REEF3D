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

#ifndef REINI_FSM_H_
#define REINI_FSM_H_

class lexer;
class fdm;
class ghostcell;
class field;

// The level set of the start made a signed distance: fast sweeping with closest-point propagation
// from the cells next to the interface (reini_fsm.cpp), called by reini_RK3 at count 0 before a few
// iterations of the reinitialisation PDE.  Returns the number of outer iterations of the sweeping.
int reini_fsm(lexer *p, fdm *a, ghostcell *pgc, field &f);

#endif
