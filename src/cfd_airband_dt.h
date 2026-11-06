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

#ifndef CFD_AIRBAND_DT_H_
#define CFD_AIRBAND_DT_H_

#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"increment.h"
#include<cmath>

// N 57 > 0 (with the level set): the time step of the convection in two parts.
// - The water and a band of N 57 cells into the air (phi on the face > -N 57 dx) keep the criterion
//   of the time step (N 47): their largest velocities go into it (uw, vw, ww).
// - The air beyond the band only needs the stability of the convection: dt <= N 58 / max(|u|/dx +
//   |v|/dy + |w|/dz) over its faces (dta; the sum of the three directions bounds the unsplit Godunov
//   convection of N 40 5/6 and the explicit schemes alike).
// Isolated velocity spikes in the air away from the free surface then no longer set the time step of
// the whole domain at the accuracy criterion of the free surface.
struct cfd_airband
{
    double uw=0.0, vw=0.0, ww=0.0;  // largest velocities in the water and the band
    double ua=0.0;                  // largest |u|, |v|, |w| in the air beyond the band (report)
    double rate=0.0;                // largest |u|/dx + |v|/dy + |w|/dz in the air beyond the band
    double dta=1.0e10;              // the time step of the air beyond the band
};

// the loops need the indices of increment
class cfd_airband_dt : public increment
{
public:
cfd_airband compute(lexer *p, fdm *a, ghostcell *pgc)
{
    cfd_airband r;

    // per cell of the air beyond the band: the rate |u|/dx + |v|/dy + |w|/dz with the face
    // velocities of the cell
    const double band = p->N57;

    ULOOP
    {
    const double ph = 0.5*(a->phi(i,j,k) + a->phi(i+1,j,k));
    const double dxl = p->DXP[IP];
    if(ph > -band*dxl)
    r.uw = MAX(r.uw, fabs(a->u(i,j,k)));
    }

    if(p->j_dir==1)
    VLOOP
    {
    const double ph = 0.5*(a->phi(i,j,k) + a->phi(i,j+1,k));
    const double dyl = p->DYP[JP];
    if(ph > -band*dyl)
    r.vw = MAX(r.vw, fabs(a->v(i,j,k)));
    }

    WLOOP
    {
    const double ph = 0.5*(a->phi(i,j,k) + a->phi(i,j,k+1));
    const double dzl = p->DZP[KP];
    if(ph > -band*dzl)
    r.ww = MAX(r.ww, fabs(a->w(i,j,k)));
    }

    LOOP
    {
    const double dl = MIN3(p->DXN[IP], (p->j_dir==1) ? p->DYN[JP] : 1.0e10, p->DZN[KP]);

        if(a->phi(i,j,k) <= -band*dl)
        {
        const double uc = MAX(fabs(a->u(i,j,k)), fabs(a->u(i-1,j,k)));
        const double vc = (p->j_dir==1) ? MAX(fabs(a->v(i,j,k)), fabs(a->v(i,j-1,k))) : 0.0;
        const double wc = MAX(fabs(a->w(i,j,k)), fabs(a->w(i,j,k-1)));

        r.rate = MAX(r.rate, uc/p->DXN[IP] + ((p->j_dir==1) ? vc/p->DYN[JP] : 0.0) + wc/p->DZN[KP]);
        r.ua = MAX(r.ua, MAX3(uc,vc,wc));
        }
    }

    double g[5] = {r.uw, r.vw, r.ww, r.ua, r.rate};
    pgc->globalmax(g,5);
    r.uw=g[0]; r.vw=g[1]; r.ww=g[2]; r.ua=g[3]; r.rate=g[4];

    if(r.rate>1.0e-20)
    r.dta = p->N58/r.rate;

    return r;
}
};

#endif
