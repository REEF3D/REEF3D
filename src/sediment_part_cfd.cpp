/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Authors: Alexander Hanke, Hans Bihs
--------------------------------------------------------------------*/

#include"sediment_part.h"
#include"partres.h"
#include"lexer.h"
#include"ghostcell.h"
#include"bedshear.h"
#include"sediment_fdm.h"
#include"bedslope.h"
#include"bedshear_reduction.h"

void sediment_part::sediment_algorithm_cfd(lexer* p, fdm* a, ghostcell* pgc, ioflow* pflow, reinitopo* preto)
{
    double starttime=pgc->timer();

    // sediment
    fill_PQ_cfd(p,a,pgc);
    pslope->slope_cds(p,pgc,s);
    pbedshear->taubed(p,a,pgc,s);
    preduce->start(p,pgc,s);
    pgc->gcsl_start4(p,s->tau_eff,1);
    pbedshear->taucritbed(p,a,pgc,s);
    pgc->gcsl_start4(p,s->tau_crit,1);

    pst->timestep(p,pgc);
    pst->move_RK2(p,a,pgc,s,pturb);
    pst->update(p,a,pgc,s,por,d50);
    pst->print_particles(p,s);

    /// topo update
    update_cfd(p,a,pgc,pflow,preto);

    p->sedsimtime=pgc->timer()-starttime;

    ++p->sediter;
}
