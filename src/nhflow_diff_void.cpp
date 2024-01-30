/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"nhflow_diff_void.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"solver.h"

nhflow_diff_void::nhflow_diff_void(lexer* p)
{
    gcval_u=10;
	gcval_v=11;
	gcval_w=12;
    
    gcval_uh=20;
	gcval_vh=21;
	gcval_wh=22;
}

nhflow_diff_void::~nhflow_diff_void()
{
}


void nhflow_diff_void::diff_u(lexer *p, fdm_nhf *d, ghostcell *pgc, solver *psolv, double *UHdiff, double *UHin, double *UH, double *VH, double *WH, double alpha)
{
    LOOP
    UHdiff[IJK] = UHin[IJK];
    
    pgc->start4V(p,UHdiff,gcval_uh);
}

void nhflow_diff_void::diff_v(lexer *p, fdm_nhf *d, ghostcell *pgc, solver *psolv, double *VHdiff, double *VHin, double *UH, double *VH, double *WH, double alpha)
{
    LOOP
    VHdiff[IJK] = VHin[IJK];
    
    pgc->start4V(p,VHdiff,gcval_vh);
}

void nhflow_diff_void::diff_w(lexer *p, fdm_nhf *d, ghostcell *pgc, solver *psolv, double *WHdiff, double *WHin, double *UH, double *VH, double *WH, double alpha)
{
    LOOP
    WHdiff[IJK] = WHin[IJK];
    
    pgc->start4V(p,WHdiff,gcval_wh);
}

void nhflow_diff_void::diff_scalar(lexer *p, fdm_nhf *d, ghostcell *pgc, solver *psolv, double *UHdiff, double *UH_in, double *UH, double *VH, double *WH, double alpha)
{
    
}