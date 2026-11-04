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

#include "ghostcell.h"
#include "lexer.h"
#include "part.h"
#include "slice.h"

void part::xchange_fillback_flag(lexer *p, ghostcell *pgc, slice &bedch, int mode)
{
    index_empty = index_empty0;

    // fill recv into F
    n = 0;
    for(int qn = 0; qn < 6; ++qn)
        for(q = 0; q < recvnum[qn]; ++q)
        {
            --index_empty;
            n = Empty[index_empty];

            // flag
            Flag[n] = ACTIVE;

            // parallel periodic boundaries: the parcel crossed the periodic side, move it by the period
            // (old and new position together, see CPM_periodic.cpp)
            if(p->periodic1 == 2)
            {
                double Lx = p->global_xmax - p->global_xmin;
                double xt = mode == 1 ? XRK1[n] : X[n];
                double sh = xt >= p->global_xmax ? -Lx : (xt < p->global_xmin ? Lx : 0.0);
                X[n] += sh;
                XRK1[n] += sh;
            }

            bedch(i, j) += ParcelFactor * Wt[n];
        }
}
