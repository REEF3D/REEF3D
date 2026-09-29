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

#include"sflow_hydrostatic.h"
#include"lexer.h"
#include"fdm2D.h" 
#include"ghostcell.h"
#include"solver2D.h"
#include"ioflow.h"
#include"patchBC_interface.h"

sflow_hydrostatic::sflow_hydrostatic(lexer* p, fdm2D *b, patchBC_interface *ppBC)
{
    pBC = ppBC;
}

sflow_hydrostatic::~sflow_hydrostatic()
{
}

void sflow_hydrostatic::start(lexer *p, fdm2D *b, ghostcell *pgc, solver2D *psolv, ioflow *pflow, 
                              slice &UH, slice &VH, slice &WH, slice &WL, slice &Un, slice &Vn, double alpha)
{
}

void sflow_hydrostatic::upgrad(lexer*p, fdm2D* b, slice &eta)
{
    SLICELOOP4
    WETDRY
    b->F(i,j) += fabs(p->W22)*eta(i,j)*(b->dfx(i,j) - b->dfx(i-1,j))/p->DXN[IP];
}

void sflow_hydrostatic::vpgrad(lexer*p, fdm2D* b, slice &eta)
{
    SLICELOOP4
    WETDRY
    b->G(i,j) += fabs(p->W22)*eta(i,j)*(b->dfy(i,j) - b->dfy(i,j-1))/p->DYN[JP];
}

void sflow_hydrostatic::wpgrad(lexer*p, fdm2D* b, slice &eta)
{
}
