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
Authors: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#include"sediment_part.h"
#include"CPM.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"vrans.h"
#include"ioflow.h"
#include"reinitopo.h"
#include"sediment_fdm.h"

// the bed seen by the fluid follows the parcels:
// topo level set from the solid volume fraction, S 10 1: bed as solid boundary,
// S 10 2: VRANS porosity 1-theta (one-way) or drag of the parcels (two-way, Q 50 1)
void sediment_part::update_cfd(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow, reinitopo *preto)
{
    pst->topo_update(p,a,pgc,s);
    preto->start(p,a,pgc,a->topo);
    pgc->start4a(p,a->topo,150);
    
    pst->bedzh_update(p,a,pgc,s);
    
    if(p->S10==1)
    {
        if(p->D22==1)
        pgc->solid_forcing_flag_update(p,a);
        
        pgc->gcdf_update(p,a);
    }
    
    pst->update(p,a,pgc,s,por,d50);
    
    // S 10 2 one-way: the fluid sees the bed as VRANS porosity
    // S 10 2 two-way (Q 50 1): the fluid sees the parcels through their drag, no porosity
    if(p->S10==2 && p->Q50!=1)
    pvrans->sedpart_update(p,a,pgc,por,d50);
    
    pflow->gcio_update(p,a,pgc);
    
    pgc->start1(p,a->u,10);
    pgc->start2(p,a->v,11);
    pgc->start3(p,a->w,12);
}
