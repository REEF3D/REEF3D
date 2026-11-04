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

#include "fdm_nhf.h"
#include "ghostcell.h"
#include "lexer.h"
#include "nhflow_forcing.h"
#include "vrans.h"

void nhflow_forcing::solid_forcing(lexer *p, fdm_nhf *d, ghostcell *pgc,
                                   double alpha, double *U, double *V, double *W, slice &WL)
{
    double du, dv, dw, udotn;
    double beta;  // 0.0 = free slip (tangential preserved), 1.0 = no slip

    if(p->A517 == 0)
        beta = 1.0;

    if(p->A517 == 1)
        beta = 0.0;

    pgc->start5V(p, d->FHB, 1);

    uf = vf = wf = 0.0;

    pgc->start5V(p, d->FHB, 1);

    uf = vf = wf = 0.0;

    // The forcing is added as  U += alpha dt CPORNH F  (nhflow_forcing, reforcing), like all momentum tendencies.
    // The direct forcing is a constraint, not a tendency: F is divided by CPORNH here, so that the solid velocity
    // is fully imposed also where a solid lies inside a porous structure (e.g. a pile with a porous cover, n < 1).
    // Without porous media CPORNH = 1 and nothing changes.

    // -- full no-slip direct forcing
    if(p->A517 == 0)
        LOOP
        {
            H = Hsolidface(p, d, 0, 0, 0);

            FX[IJK] += H * (uf - U[IJK]) / (alpha * p->dt * CPORNH);
            FY[IJK] += H * (vf - V[IJK]) / (alpha * p->dt * CPORNH);
            FZ[IJK] += H * (wf - W[IJK]) / (alpha * p->dt * CPORNH);

            d->FHB[IJK] = MIN(d->FHB[IJK] + H, 1.0);
        }

    // -- no-penetration in the band, tangential left to the wall model
    if(p->A517 == 1)
        LOOP
        {
            H = Hsolidface(p, d, 0, 0, 0);

            d->FHB[IJK] = MIN(d->FHB[IJK] + H, 1.0);

            du = uf - U[IJK];
            dv = vf - V[IJK];
            dw = wf - W[IJK];

            if(d->SOLID[IJK] <= 0.0)
            {
                // interior of the body: enforce the full rigid-body velocity
                FX[IJK] += H * du / (alpha * p->dt * CPORNH);
                FY[IJK] += H * dv / (alpha * p->dt * CPORNH);
                FZ[IJK] += H * dw / (alpha * p->dt * CPORNH);
            }

            // -- no-penetration in the band, tangential left to the wall model
            if(p->A517 == 1)
                LOOP
                {
                    H = Hsolidface(p, d, 0, 0, 0);

                    d->FHB[IJK] = MIN(d->FHB[IJK] + H, 1.0);

                    du = uf - U[IJK];
                    dv = vf - V[IJK];
                    dw = wf - W[IJK];

                    if(d->SOLID[IJK] <= 0.0)
                    {
                        nx /= norm;
                        ny /= norm;
                        nz /= norm;

                        udotn = nx * du + ny * dv + nz * dw;

                        FX[IJK] += H * (beta * du + (1.0 - beta) * udotn * nx) / (alpha * p->dt * CPORNH);
                        FY[IJK] += H * (beta * dv + (1.0 - beta) * udotn * ny) / (alpha * p->dt * CPORNH);
                        FZ[IJK] += H * (beta * dw + (1.0 - beta) * udotn * nz) / (alpha * p->dt * CPORNH);
                    }

                    else
                    {
                        // degenerate gradient (flat level set): fall back to full forcing
                        FX[IJK] += H * du / (alpha * p->dt * CPORNH);
                        FY[IJK] += H * dv / (alpha * p->dt * CPORNH);
                        FZ[IJK] += H * dw / (alpha * p->dt * CPORNH);
                    }
                }

            pgc->start5V(p, d->FHB, 50);
        }