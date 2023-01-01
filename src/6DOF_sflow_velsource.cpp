/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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
Author: Tobias Martin
--------------------------------------------------------------------*/

#include"6DOF_sflow.h"
#include"lexer.h"
#include"fdm.h"
#include"fdm2D.h"
#include"ghostcell.h"

void sixdof_sflow::isource(lexer *p, fdm *a, ghostcell *pgc)
{
}

void sixdof_sflow::jsource(lexer *p, fdm *a, ghostcell *pgc)
{
}

void sixdof_sflow::ksource(lexer *p, fdm *a, ghostcell *pgc)
{
}

void sixdof_sflow::isource2D(lexer *p, fdm2D *b, ghostcell *pgc)
{
	SLICELOOP1
    {
        b->F(i,j) += 1.0/p->W1*(press(i+1,j) - press(i,j))/p->DXP[IP];
    }
}

void sixdof_sflow::jsource2D(lexer *p, fdm2D *b, ghostcell *pgc)
{
	SLICELOOP2
    {
        b->G(i,j) += 1.0/p->W1*(press(i,j+1) - press(i,j))/p->DYP[JP];
    }
}
