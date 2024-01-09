/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"ioflow_f.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"

void ioflow_f::ini(lexer *p, fdm* a, ghostcell* pgc)
{
    gcio_update(p,a,pgc);
}

void ioflow_f::ini_fnpf(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
}

void ioflow_f::ini_ptf(lexer *p, fdm* a, ghostcell* pgc)
{
}

void ioflow_f::ini2D(lexer *p, fdm2D* b, ghostcell* pgc)
{
    discharge2D(p,b,pgc);
    inflow2D(p,b,pgc,b->P,b->Q,b->bed,b->eta);
}

void ioflow_f::full_initialize2D(lexer *p, fdm2D *b, ghostcell *pgc)
{
}

