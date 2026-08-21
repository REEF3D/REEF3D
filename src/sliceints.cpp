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

#include "sliceint.h"
#include "sliceint1.h"
#include "sliceint2.h"
#include "sliceint4.h"
#include "sliceint5.h"

// Out-of-line destructors: each is its class's key function, so the vtable
// is emitted once, here, instead of in every translation unit that includes
// the header (-Wweak-vtables).

sliceint::~sliceint() = default;
sliceint1::~sliceint1() = default;
sliceint2::~sliceint2() = default;
sliceint4::~sliceint4() = default;
sliceint5::~sliceint5() = default;
