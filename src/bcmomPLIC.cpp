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

#include "bcmom.h"
#include "lexer.h"
#include "fdm.h"
#include "ghostcell.h"
#include "turbulence.h"
#include "VOF_PLIC.h"
#include "bc_noflux.h"

void bcmom::bcmomPLIC_start(fdm* a, lexer* p, ghostcell *pgc, turbulence *pturb, VOF_PLIC *pplic, field& b, int gcval)
{
    wall_laws(p,a,pturb,b,gcval);

    pplic->surface_tension2D(p,a,pgc,gcval);
}
