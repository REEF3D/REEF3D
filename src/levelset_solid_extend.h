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

#ifndef LEVELSET_SOLID_EXTEND_H_
#define LEVELSET_SOLID_EXTEND_H_

class lexer;
class fdm;
class ghostcell;
class field;

// The level set inside the solids of the direct forcing (S entities of the grid file, G 505 0) extrapolated
// from the fluid: constant across, linear along z (X 48 2, the default; levelset_solid_extend.cpp)
void levelset_solid_extend(lexer *p, fdm *a, ghostcell *pgc, field &f, int gcval);

#endif
