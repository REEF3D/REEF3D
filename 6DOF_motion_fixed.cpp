/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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
--------------------------------------------------------------------*/

#include"6DOF_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void sixdof_f::motion_fixed(lexer *p, fdm *a, ghostcell *pgc)
{
	if(p->X210==1)
	{
	Uext = p->X210_u;
	Vext = p->X210_v;
	Wext = p->X210_w;
	}
	
	if(p->X211==1)
	{
	Pext = p->X211_p;
	Qext = p->X211_q;
	Rext = p->X211_r;
	}
}