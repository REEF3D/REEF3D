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

#include"fdm_seastate.h"
#include"lexer.h"

fdm_seastate::fdm_seastate(lexer *p) : bed(p),depth(p),eta(p),U(p),V(p),
                                      Hs(p),Tm01(p),Tm10(p),Tp(p),dir(p),spread(p),
                                      ddx(p),ddy(p),dUdx(p),dUdy(p),dVdx(p),dVdy(p),dddt(p),depth_n(p),
                                      wet(p),wet0(p),refr(p),nodeval(p)
{
}
