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

#include "fnpf_scour.h"
#include "fdm_fnpf.h"
#include "fnpf_sigma_col.h"
#include "ghostcell.h"
#include "lexer.h"
#include <cmath>

fnpf_scour::fnpf_scour(lexer *p, fdm_fnpf *c, ghostcell *pgc) : scour_core(p, pgc, "FNPF", &c->depth, false)
{
}

void fnpf_scour::sample(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    bedvel(p, c);
    scour_core::sample(p, pgc);
}

void fnpf_scour::bedvel(lexer *p, fdm_fnpf *c)
{
    // the horizontal velocity at the bed node (k = 0): u = Fi_x + sigx Fi_sig, with Fi_sig from the
    // bed boundary condition of the solved potential (as the ghost cells of fnpf_bed_update::bedbc_sig)
    const double *const Fi = c->Fi;
    const double ydir = p->y_dir;

    k = 0;
    SLICELOOP4
    {
        if(p->wet[IJ] == 0 || p->flag7[FIJK] <= 0)
        {
            ub(i, j) = vb(i, j) = 0.0;
            continue;
        }

        const double Fx = (Fi[FIp1JK] - Fi[FIm1JK]) / (p->DXP[IP] + p->DXP[IM1]);
        const double Fy = ydir * (Fi[FIJp1K] - Fi[FIJm1K]) / (p->DYP[JP] + p->DYP[JM1]);

        double sgx, sgy;

        if(p->sigx != nullptr)
        {
            sgx = p->sigx[FIJK];
            sgy = p->sigy[FIJK];
        }
        else
        {
            fnpf_sigma_col col;
            col.set(p, c, i, j);
            sgx = col.sx(p->ZN[KP]);
            sgy = col.sy(p->ZN[KP]);
        }

        const double bcval = c->Bx(i, j) * Fx + c->By(i, j) * Fy;
        double denom = p->sigz[IJ] + c->Bx(i, j) * sgx + c->By(i, j) * sgy;

        if(fabs(denom) < 1.0e-20)
            denom = 1.0e-20;

        const double Fs = -bcval / denom;

        ub(i, j) = Fx + sgx * Fs;
        vb(i, j) = ydir * (Fy + sgy * Fs);
    }
}
