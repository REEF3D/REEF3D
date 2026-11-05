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

#ifndef PATCHBC_CODES_H_
#define PATCHBC_CODES_H_

// Boundary codes of the patch faces (gcb[n][4], gcbsl[n][4]), set by patchBC_ini on the wall faces
// (21/22) inside a patch geometry B 440-442. Each patch is exactly one of two well-posed kinds:
//
//   inlet  (B 411, 414, 415, 421): velocity Dirichlet (set by patchBC), pressure Neumann,
//                                  pressure correction Neumann
//   outlet (B 412, 413, 418, 422 or nothing): velocity zero gradient, pressure Dirichlet
//                                  (hydrostatic, set by patchBC), pressure correction Dirichlet 0
//
// _FSF: the level set / water level in the ghost cells is set by patchBC (B 413, B 422),
// otherwise it is extrapolated by the ghostcell routines.
// Use the predicates below, never the numbers.

enum {
    PATCH_INLET = 61,
    PATCH_OUTLET = 62,
    PATCH_INLET_FSF = 63,
    PATCH_OUTLET_FSF = 64
};

inline bool patch_bc(int bc)
{
    return bc>=PATCH_INLET && bc<=PATCH_OUTLET_FSF;
}

inline bool patch_inlet(int bc)
{
    return bc==PATCH_INLET || bc==PATCH_INLET_FSF;
}

inline bool patch_outlet(int bc)
{
    return bc==PATCH_OUTLET || bc==PATCH_OUTLET_FSF;
}

inline bool patch_fsf(int bc)
{
    return bc==PATCH_INLET_FSF || bc==PATCH_OUTLET_FSF;
}

#endif
