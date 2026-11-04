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

#ifndef NHFLOW_BREAKING_PROFILE_H_
#define NHFLOW_BREAKING_PROFILE_H_

#include"lexer.h"
#include"fdm_nhf.h"
#include<cmath>

// breaking viscosity vb of column (ii,jj) at the relative height sig (0 bed, 1 surface).
// A 550 4: vb cosh(k sig h)/cosh(k h), the profile of the orbital velocity with the wavenumber k of the breaking
// crest (d->vbk, 2 pi / trough-to-trough distance) and the water depth h: uniform in shallow water, exp(-k (h - z))
// in deep water. Otherwise vb.
inline double nhflow_vbz(lexer *p, fdm_nhf *d, int ii, int jj, double sig)
{
    const double v = d->vb(ii,jj);

    if(p->A550!=4 || v<=0.0)
    return v;

    const double kh = d->vbk(ii,jj)*MAX(d->WL(ii,jj),0.0);

    if(kh<=0.0)
    return v;

    // cosh(kh s)/cosh(kh) without overflow
    return v*(exp(kh*(sig-1.0)) + exp(-kh*(sig+1.0)))/(1.0 + exp(-2.0*kh));
}

#endif
