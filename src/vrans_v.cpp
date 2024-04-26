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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"vrans_v.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

vrans_v::vrans_v(lexer *p, ghostcell *pgc) 
{
}

vrans_v::~vrans_v()
{
}

void vrans_v::veltimesave(lexer *p, fdm *a, ghostcell *pgc)
{
}

void vrans_v::initialize(lexer *p, fdm *a, ghostcell *pgc)
{
}

void vrans_v::sed_update(lexer *p, fdm *a, ghostcell *pgc)
{
}

void vrans_v::u_source(lexer *p, fdm *a)
{
}

void vrans_v::v_source(lexer *p, fdm *a)
{
}

void vrans_v::w_source(lexer *p, fdm *a)
{
}

void vrans_v::kw_source(lexer *p, fdm *a, field &kin)
{
}

void vrans_v::ke_source(lexer *p, fdm *a, field &kin)
{
}

void vrans_v::eps_source(lexer *p, fdm *a, field &kin, field &eps)
{
}

void vrans_v::omega_source(lexer *p, fdm *a, field &kin, field &eps)
{
}

void vrans_v::eddyv_func(lexer *p, fdm *a)
{
}
