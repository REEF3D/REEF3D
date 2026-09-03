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

#ifndef FIELDINT7_H_
#define FIELDINT7_H_

#include "fieldint.h"

// Integer counterpart of field7: sigma-grid vertical-node layout, stride
// p->kmaxF, addressed with the FIJK family. See field7.h for the slack plane.
class fieldint7 final : public fieldint
{
public:
    fieldint7(lexer *p) : fieldint(p, p->kmaxF,
                                   static_cast<std::size_t>(p->imax)*static_cast<std::size_t>(p->jmax)) {}
    virtual ~fieldint7() = default;
};

#endif
