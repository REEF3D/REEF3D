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

#ifndef CFD_CKPT_H_
#define CFD_CKPT_H_

#include<string>

class lexer;
class reefckpt;

// CFD in checkpoints (P 48, I 45): the fields and slices of the level-0 grid (cfd_ckpt.cpp); NHFLOW and SFLOW
// with sediment the same for their grid (the slices of the sediment objects), under another tag
void cfd_ckpt_write(lexer*, reefckpt&, const std::string &tag="cfd");
bool cfd_ckpt_read(lexer*, reefckpt&, std::string&, const std::string &tag="cfd");

#endif
