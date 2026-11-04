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

#ifndef IOFLOW_GCIO_H_
#define IOFLOW_GCIO_H_

class lexer;
class patchBC_interface;

// shared part of gcio_update (ioflow_f, iowave, ioflow_gravity; CFD and NHFLOW)
//
// inflow  boundary cells: codes 1 (inflow) and 6 (wave generation)
// outflow boundary cells: codes 2 (outflow), 7 (numerical beach), 8 (outflow 2)

// rebuild p->gcin / p->gcout from the boundary cells with flag[IJK]>0
// (and wet[IJ]==1 for the inflow cells if wet is given)
void ioflow_gcio_lists(lexer *p, const int *flag, const int *wet);

// reset p->IO and mark the ghost cells: 1 inflow, 2 outflow (cells with flag[IJK]>0), 3 patch faces
void ioflow_gcio_marks(lexer *p, const int *flag, patchBC_interface *pBC);

#endif
