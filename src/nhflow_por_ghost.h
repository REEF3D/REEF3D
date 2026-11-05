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

#ifndef NHFLOW_POR_GHOST_H_
#define NHFLOW_POR_GHOST_H_

#include"lexer.h"
#include"fdm_nhf.h"

// NHFLOW porosity POR in the cells below the bed and above the surface of EVERY column of the array:
// own columns, partition halo, ghost columns beyond the domain and the 2D y-ghost columns.
// start5Vfull fills the own columns only. The Poisson rows of the bottom layer read POR below the bed
// of the neighbour column (cN, cS, cW, cE of nhflow_poisson.cpp); the ghost columns beyond the domain
// and the 2D y-ghost columns were never written and stayed 0, so CPORNHval(0) = 0/0 = NaN with B 260 0
// (Poisson residual NaN in the first step) and 0 (coefficient 0 instead of 1) otherwise.
// Called after every update of POR (initialisation, VRANS update), after the halo exchange.
inline void nhflow_por_ghost(lexer *p, fdm_nhf *d)
{
    for(int ii=p->imin; ii<p->imin+p->imax; ++ii)
    for(int jj=p->jmin; jj<p->jmin+p->jmax; ++jj)
    {
        const int c0 = (ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax - p->kmin;

        // fluid layers that were never written (ghost and corner columns): open flow
        for(int kk=0; kk<p->knoz; ++kk)
        if(!(d->POR[c0+kk]>0.0))
        d->POR[c0+kk] = 1.0;

        const double nb = d->POR[c0];
        const double nt = d->POR[c0+p->knoz-1];

        for(int kk=p->kmin; kk<0; ++kk)
        d->POR[c0+kk] = nb;

        for(int kk=p->knoz; kk<p->kmin+p->kmax; ++kk)
        d->POR[c0+kk] = nt;
    }
}

#endif
