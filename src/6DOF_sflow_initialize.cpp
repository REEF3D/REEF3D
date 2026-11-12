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
Authors: Hans Bihs, Tobias Martin
--------------------------------------------------------------------*/

#include"6DOF_sflow.h"
#include"lexer.h"
#include"fdm.h"
#include"fdm2D.h"
#include"ghostcell.h"
   
void sixdof_sflow::ini(lexer *p, ghostcell *pgc)
{
}

void sixdof_sflow::initialize(lexer *p, fdm2D *b, ghostcell *pgc)
{
    for (int nb = 0; nb < number6DOF; nb++)
    fb_obj[nb]->initialize_shipwave(p, pgc,b->eta,b->WL);
    
    // X 402 1: the pressure of the initial position is active from the first step
    if(p->X10==3 && p->X402==1)
    {
    update_press(p,pgc,b->eta);
    pressure_head(p,b,pgc);
    
    // no draft ramp (X 207 0): start from the equilibrium depression of the applied pressure,
    // eta = pbh, for all hull types (X 400 2, 3: the analytic pressure, not the STL draft)
    if(p->X207==0)
    SLICELOOP4
    {
    b->eta(i,j) = b->pbh(i,j);
    b->WL(i,j) = MAX(b->eta(i,j) + p->wd - p->bed[IJ],0.0);
    }
    }
}

void sixdof_sflow::initialize(lexer *p, fdm *a, ghostcell *pgc)
{
}

void sixdof_sflow::initialize(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
}


void sixdof_sflow::start_cfd(lexer* p, fdm* a, ghostcell* pgc, int iter, 
                                field &uvel, field &vvel, field &wvel, field &fx, field &fy, field &fz, bool finalize)
{
}

void sixdof_sflow::start_nhflow(lexer* p, fdm_nhf* d, ghostcell* pgc, int iter, 
                                        double *U, double *V, double *W, double *FX, double *FY, double *FZ, slice &WL, slice &fe, bool finalize)
{
}