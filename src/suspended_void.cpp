/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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

#include"suspended_void.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"diffusion.h"
#include"ioflow.h"
#include"turbulence.h"
#include"solver.h"

suspended_void::suspended_void()
{
}

suspended_void::~suspended_void()
{
}

void suspended_void::start(fdm* a, lexer* p, convection* pconvec, diffusion* pdiff, solver* psolv, ghostcell* pgc, ioflow* pflow)
{
}

void suspended_void::ctimesave(lexer *p, fdm* a)
{
}
