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
for more details->

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Authors: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#include"sediment_part.h"
#include"CPM.h"
#include"lexer.h"
#include"ghostcell.h"
#include"fdm.h"
#include"bedshear.h"
#include"sediment_fdm.h"

void sediment_part::ini_cfd(lexer *p, fdm *a, ghostcell *pgc)
{
    double h;
    ILOOP
    JLOOP
    {
        h = p->ZN[0+marge];
        
        KLOOP
        PBASECHECK
        if(a->topo(i,j,k-1)<0.0 && a->topo(i,j,k)>0.0)
            h = -(a->topo(i,j,k-1)*p->DZP[KP])/(a->topo(i,j,k)-a->topo(i,j,k-1)) + p->pos_z()-p->DZP[KP];

        s->bedzh(i,j)=h;
        s->bedzh0(i,j)=h;
    }

    pgc->gcsl_start4(p,s->bedzh0,50);

    pst->seed_particles(p,a,pgc,s);

    fill_PQ_cfd(p,a,pgc);

    SLICELOOP4
        s->bedk(i,j)=0;

    SLICELOOP4
    {
        KLOOP
        PBASECHECK
        if(a->topo(i,j,k)<0.0 && a->topo(i,j,k+1)>=0.0)
        s->bedk(i,j)=k+1;

        s->reduce(i,j)=0.3;
    }

    pbedshear->taubed(p,a,pgc,s);
    pgc->gcsl_start4(p,s->tau_eff,1);
    pbedshear->taucritbed(p,a,pgc,s);
    pgc->gcsl_start4(p,s->tau_crit,1);

    pst->ini_fields(p,a,pgc,s);
    pst->print_particles(p,s);
    
    // S 10 2 two-way (Q 50 1): the bed is no solid for the fluid, the fluid flows through
    // the parcels and feels their drag; topo stays as the bed level for output
    if(p->S10==2 && p->Q50==1 && p->Q11==2)
    {
        p->toporead=0;
        p->topoforcing=0;
        
        int G5 = p->G5;
        p->G5 = 1;
        pgc->gcdf_update(p,a);
        p->G5 = G5;
    }

    pgc->gcdf_update(p,a);
}
