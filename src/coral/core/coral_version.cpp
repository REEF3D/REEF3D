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

#include"coral/core/coral_version.h"

#define CORAL_NAME "REEF3D::CORAL"
#define CORAL_LONG_NAME "Coupled Offshore Rotor Analysis and Loads"
#define CORAL_VERSION "0.1.0"

namespace coral
{
    const char *name() {return CORAL_NAME;}
    const char *long_name() {return CORAL_LONG_NAME;}
    const char *version() {return CORAL_VERSION;}
    const char *banner() {return CORAL_NAME " " CORAL_VERSION " (" CORAL_LONG_NAME ")";}
}
