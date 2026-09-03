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
Author: Alexander Hanke
--------------------------------------------------------------------*/

#ifndef FIELD7_H_
#define FIELD7_H_

#include "field.h"

// Sigma-grid vertical-node field: n*m*(o+1) planes, vertical stride p->kmaxF.
//
// Same (i,j,k) accessor as field1-4 - only the vertical extent differs. The
// vertical index runs 0..knoz+1 (FKLOOP) rather than 0..knoz (KLOOP), i.e. one
// plane per cell plus a closing top plane, which is what the sigma-grid solvers
// (FNPF, NHFLOW) need. The flat-array equivalent is the FIJK family in
// iterators3D.h; field7 exists so those call sites can drop the flat index and
// address by (i,j,k) like every other field.
class field7 : public field
{
public:
    // The slack plane mirrors the imax*jmax*(kmax+2) allocation the sigma grid
    // makes against a kmaxF = kmax+1 stride. FIJKp3 at the top node of the last
    // column lands on the very last in-stride slot, so the spare plane is the
    // only headroom the forward-stencil macros have. Drop it only after
    // auditing FIJKp3/p4.
    field7(lexer *p) : field(p, p->kmaxF,
                             static_cast<std::size_t>(p->imax)*static_cast<std::size_t>(p->jmax)) {}
    virtual ~field7() = default;
};

#endif
