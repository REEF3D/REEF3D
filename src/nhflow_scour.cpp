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

#include "nhflow_scour.h"
#include "fdm_nhf.h"
#include "ghostcell.h"
#include "lexer.h"

nhflow_scour::nhflow_scour(lexer *p, fdm_nhf *d, ghostcell *pgc) : scour_core(p, pgc, "NHFLOW", &d->depth, true)
{
    if(p->geo_curv == 1 && (p->P196 == 1 || p->P198 > 0) && p->mpirank == 0)
        cout << "NHFLOW scour: the current file P 196 and the gauges P 198 locate points on the Cartesian grid, "
             << "not on a curvilinear grid" << endl;
}

void nhflow_scour::sample(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    SLICELOOP4
    {
        ub(i, j) = vb(i, j) = ucx(i, j) = ucy(i, j) = 0.0;

        if(p->wet[IJ] == 0)
            continue;

        // near bed: the lowest layer
        k = 0;
        ub(i, j) = d->U[IJK];
        vb(i, j) = d->V[IJK];

        // depth average
        double us = 0.0, vs = 0.0, ws = 0.0;

        for(k = 0; k < p->knoz; ++k)
        {
            const double w = p->DZN[KP];

            us += w * d->U[IJK];
            vs += w * d->V[IJK];
            ws += w;
        }

        if(ws > 1.0e-12)
        {
            ucx(i, j) = us / ws;
            ucy(i, j) = vs / ws;
        }
    }

    scour_core::sample(p, pgc);
}
