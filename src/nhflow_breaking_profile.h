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
#include"increment.h"
#include<cmath>

// breaking viscosity vb of column (ii,jj) at the relative height sig (0 bed, 1 surface).
// A 550 4, 5: vb cosh(k sig h)/cosh(k h), the profile of the orbital velocity with the wavenumber k of the breaking
// crest (d->vbk, 2 pi / trough-to-trough distance) and the water depth h: uniform in shallow water, exp(-k (h - z))
// in deep water. Otherwise vb.
inline double nhflow_vbz(lexer *p, fdm_nhf *d, int ii, int jj, double sig)
{
    const double v = d->vb(ii,jj);

    if((p->A550!=4 && p->A550!=5) || v<=0.0)
    return v;

    const double kh = d->vbk(ii,jj)*MAX(d->WL(ii,jj),0.0);

    if(kh<=0.0)
    return v;

    // cosh(kh s)/cosh(kh) without overflow
    return v*(exp(kh*(sig-1.0)) + exp(-kh*(sig+1.0)))/(1.0 + exp(-2.0*kh));
}

// A 556: share of the breaking viscosity on the faces between the sigma layers (vertical stress), uniform over the
// depth. The Kennedy et al. (2000) viscosity is the one of a depth-integrated model; in full on the vertical faces
// it makes the undertow depth-uniform and weak and smooths the wave bottom boundary layer (Cox et al. 1996: undertow
// -6.4 cm/s at L4, measured -11.6; with A 556 0.1: -8.1 and sheared). With no vertical share (and with a share that
// goes to 0 at the bed) the bore front grows grid-scale noise: the mixed sigma terms of the horizontal viscosity,
// 2 sigx d2u/dx dsigma, are explicit and rely on the implicit vertical part for their stability.
inline double nhflow_vbv(lexer *p)
{
    return MAX(0.0, MIN(1.0, p->A556));
}

// A 556 < 1: the horizontal part of vb acts along the sigma layers (Kennedy et al. 2000 diffuse along the layer of a
// depth-integrated model): on the horizontal faces only, without the sigx, sigy parts of the sigma faces, the sigxx
// terms and the explicit mixed terms 2 sigx d2u/dx dsigma. Those mixed terms are explicit and need an implicit
// vertical part of the same viscosity for their stability; with the full vb in them A 556 0 grew grid-scale noise at
// the bore front.
inline bool nhflow_vbsig(lexer *p)
{
    return nhflow_vbv(p)<1.0;
}

// the breaking viscosity of layer kk of column (ii,jj) as the momentum sees it: A 550 1 in the top layer, A 550 3/4/5
// in all layers (A 550 4/5 with the depth profile of nhflow_vbz)
inline double nhflow_vb_layer(lexer *p, fdm_nhf *d, int ii, int jj, int kk)
{
    if(p->A550>=3)
    return nhflow_vbz(p,d,ii,jj,p->ZP[kk+increment::marge]);

    return (p->A550>=1 && kk==p->knoz-1) ? d->vb(ii,jj) : 0.0;
}

#endif
