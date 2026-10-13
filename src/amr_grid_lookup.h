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


#ifndef AMR_GRID_LOOKUP_H_
#define AMR_GRID_LOOKUP_H_

class lexer;
class fdm;

// The finest grid of a mesh refinement at a point (rank local): the loads of a body (sixdof_obj_cfd
// forces_stl) take the pressure, density, viscosity and velocities there.  false: level 0.
class amr_grid_lookup
{
public:
    virtual ~amr_grid_lookup() {}
    virtual bool finest_grid(double x, double y, double z, lexer *&q, fdm *&a) const = 0;
};

#endif
