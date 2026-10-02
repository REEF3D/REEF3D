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

/*--------------------------------------------------------------------
DEM interface: standalone discrete element module for rigid particles
of arbitrary shape, coupled to REEF3D::CFD and REEF3D::NHFLOW.
--------------------------------------------------------------------*/

#ifndef DEM_H_
#define DEM_H_

class lexer;
class fdm;
class fdm_nhf;
class ghostcell;
class field;
class slice;

class dem
{
public:
    virtual ~dem() = default;

    // advance the particles over one fluid time step, called before the momentum step
    virtual void start_cfd(lexer*, fdm*, ghostcell*)=0;
    virtual void start_nhflow(lexer*, fdm_nhf*, ghostcell*)=0;

    // fluid forcing inside the momentum RK stages (resolved: direct forcing, unresolved: momentum source)
    virtual void forcing_cfd(lexer*, fdm*, ghostcell*, int, double, field&, field&, field&, bool)=0;
    virtual void forcing_nhflow(lexer*, fdm_nhf*, ghostcell*, int, double, double*, double*, double*, slice&, bool, bool)=0;
};

#endif
