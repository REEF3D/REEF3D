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

#include"mast/core/mast_version.h"

#define MAST_NAME "REEF3D::MAST"
#define MAST_LONG_NAME "Multibody Aero-hydro-Servo-elastic Toolkit"
#define MAST_VERSION "0.1.0"

namespace mast
{
    const char *name() {return MAST_NAME;}
    const char *long_name() {return MAST_LONG_NAME;}
    const char *version() {return MAST_VERSION;}
    const char *banner() {return MAST_NAME " " MAST_VERSION " (" MAST_LONG_NAME ")";}
}
