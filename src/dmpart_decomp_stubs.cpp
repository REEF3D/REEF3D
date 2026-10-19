/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

This file is part of REEF3D (port of the DIVEMesh decomposition, dmpart).

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

#include "dmpart_decomp.h"

// the output files of the DIVEMesh decomposition (DIVEMesh_Decomp) are not written in REEF3D
namespace dmpart
{
void decomp::print_partition(lexer*, dive*) {}
void decomp::partition_planes(lexer*, dive*) {}
void decomp::decomp_vtp(lexer*, dive*) {}
}
