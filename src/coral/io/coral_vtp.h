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
// REEF3D::CORAL: VTK PolyData (ASCII .vtp) of the structure for ParaView: beam elements and mooring
// lines as lines, the bodies as vertices (cell data "kind": 0 body, 1 beam element, 2 mooring line,
// 3 superelement element between recovered nodes).

#ifndef CORAL_VTP_H_
#define CORAL_VTP_H_

#include<string>

namespace coral
{
    class model;

    // false if the file cannot be written
    bool write_vtp(const model&, const std::string &file);
}

#endif
