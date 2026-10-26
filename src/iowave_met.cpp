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

#include"iowave.h"
#include"lexer.h"
#include"ghostcell.h"
#include"met_provider.h"

/*--------------------------------------------------------------------
Inverse barometer of the met forcing (W 200 > 0, W 210 1) on the open
boundaries of SFLOW and NHFLOW: the level the Riemann, Flather and
clamped level edges (B 520 3-5), the relaxation (wave generation) and
the beach zones target is raised by

  eta_ib = -r (pa - pref) / (rho_w g),   pref = W 209, r the ramp W 204

with the air pressure pa of the met forcing at the model time of the
step, so that a low pressure enters the domain with its static response
and is not held at still water by the boundary.
--------------------------------------------------------------------*/

double iowave::eta_ib(lexer *p, ghostcell *pgc, double x, double y)
{
    return met_provider::get(p,pgc)->ib(p,pgc,x,y);
}
