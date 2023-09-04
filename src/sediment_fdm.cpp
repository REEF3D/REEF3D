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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"sediment_fdm.h"
#include"lexer.h"

sediment_fdm::sediment_fdm(lexer *p) : P(p),Q(p),
                                       bedzh(p),bedzh0(p),bedch(p),bedsole(p),
                                       vz(p),dh(p),reduce(p),
                                       ks(p),
                                       tau_eff(p),tau_crit(p),shearvel_eff(p),shearvel_crit(p),shields_eff(p),shields_crit(p),
                                       qb(p),qbe(p),
                                       alpha(p),teta(p),gamma(p),beta(p),phi(p),
                                       active(p),
                                       bedk(p),slideflag(p),
                                       cbe(p),cb(p),cbn(p),conc(p),
                                       waterlevel(p),guard(p)
{

}

sediment_fdm::~sediment_fdm()
{
}












