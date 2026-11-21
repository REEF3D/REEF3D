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

#include"sediment_fdm.h"
#include"lexer.h"

sediment_fdm::sediment_fdm(lexer *p) : P(p),Q(p),
                                       bedzh(p),bedzh0(p),bedch(p),bedsole(p),
                                       vz(p),dh(p),reduce(p),
                                       ks(p),ks_eff(p),ro(p),
                                       tau_eff(p),tau_crit(p),shearvel_eff(p),shearvel_crit(p),shields_eff(p),shields_crit(p),
                                       qb(p),qbe(p),qbs(p),
                                       alpha(p),teta(p),gamma(p),beta(p),phi(p),
                                       active(p),
                                       bedk(p),slide_fh(p),
                                       cbe(p),cb(p),cbn(p),conc(p),rbe(p),wsx(p),rbf(p),rbz(p),ub(p),vb(p),dryd(p),
                                       waterlevel(p),guard(p),MOB(p),tau_i(p)
{
    dk = p->S20;
    pmix = nullptr;
    
    for(int n=0; n<p->imax*p->jmax; ++n)
    {
    rbe.V[n] = 1.0;
    wsx.V[n] = 0.0;
    rbf.V[n] = 1.0;
    rbz.V[n] = 0.0;
    }

}

sediment_fdm::~sediment_fdm()
{
}

// S 41 1: after the sediment step of flow step S 43 (the suspended solve of that step ran before it, without a
// reference concentration); S 41 2, 3: the start criteria of sediment_f::start_cfd / start_nhflow.
// The exchange was gated with S 43 for every S 41: with S 41 2 or 3 and S 43 at its default (1000) the water
// column picked up sediment from the first sediment step on and the bed gave none until flow step 1000.
bool sediment_fdm::started(lexer *p)
{
    if(p->S41==2)
    return p->simtime>=p->S45;
    
    if(p->S41==3)
    return p->simtime/p->wT>=p->S47 && p->count>0;
    
    return p->count>p->S43;
}
