/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Authors: Hans Bihs
--------------------------------------------------------------------*/

#include"6DOF_motionext_file.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

sixdof_motionext_file::sixdof_motionext_file(lexer *p, ghostcell *pgc)
{
    ini(p,pgc);
}

sixdof_motionext_file::~sixdof_motionext_file()
{
}

void sixdof_motionext_file::read_format_1(lexer *p, ghostcell *pgc)
{
    Uext = p->X210_u;
    Vext = p->X210_v;
    Wext = p->X210_w;

    Pext = p->X211_p;
    Qext = p->X211_q;
    Rext = p->X211_r;
}
