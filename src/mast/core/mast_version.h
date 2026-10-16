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
// REEF3D::MAST: name and version.

#ifndef MAST_VERSION_H_
#define MAST_VERSION_H_

namespace mast
{
    // "REEF3D::MAST"
    const char *name();
    // "Multibody Aero-hydro-Servo-elastic Toolkit"
    const char *long_name();
    // version of the MAST module, independent of the REEF3D git version
    const char *version();
    // one line for logs: name, version and long name
    const char *banner();
}

#endif
