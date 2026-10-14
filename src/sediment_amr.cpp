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

#include"sediment_f.h"
#include"sediment_fdm.h"
#include"sediment_exner.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"bedload.h"
#include"bedconc.h"
#include"bedshear.h"
#include"sandslide.h"
#include"topo_relax.h"
#include"bedslope.h"
#include"bedshear_reduction.h"
#include"bedload_direction.h"
#include"sediment_roughness.h"

// The sediment step of one grid of the mesh refinement or of a block, in the parts that the grids run
// before and after their shared Exner step (nhflow_amr_sed.cpp, sflow_amr_sed.cpp). The parts are those of
// sediment_algorithm_nhflow / sediment_algorithm_sflow.

sediment_exner* sediment_f::exner()
{
    return static_cast<sediment_exner*>(ptopo);
}

void sediment_f::amr_physics_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    prep_nhflow(p,d,pgc);

    pslope->slope_cds(p,pgc,s);
    preduce->start(p,pgc,s);
    pks->start(p,pgc,s,d->WL);

	pbedshear->taubed(p,d,pgc,s);
    pbedshear->taucritbed(p,d,pgc,s);

    pbed->start(p,pgc,s);
    pbeddir->start(p,pgc,s);
    pcbed->start(p,pgc,s);

	prelax->start(p,pgc,s);

    exner()->amr_qb(p,pgc,s);
}

void sediment_f::amr_physics_sflow(lexer *p, fdm2D *b, ghostcell *pgc)
{
    prep_sflow(p,b,pgc,b->P,b->Q);

    pslope->slope_cds(p,pgc,s);
    preduce->start(p,pgc,s);
    pks->start(p,pgc,s,b->hp);

	pbedshear->taubed(p,b,pgc,s);
    pbedshear->taucritbed(p,b,pgc,s);

    pbed->start(p,pgc,s);
    pbeddir->start(p,pgc,s);

	prelax->start(p,pgc,s);

    exner()->amr_qb(p,pgc,s);
}

void sediment_f::amr_finish_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, ioflow *pflow)
{
	prelax->start(p,pgc,s);

    update_nhflow(p,d,pgc,pflow);
}

void sediment_f::amr_finish_sflow(lexer *p, fdm2D *b, ghostcell *pgc, ioflow *pflow)
{
	prelax->start(p,pgc,s);

    update_sflow(p,b,pgc,pflow);
}

// ---- the sand slide in steps (sandslide.h)

void sandslide::amr_begin(lexer *p, sediment_fdm *s)
{
    int i,j;

    SEDSLICELOOP
    s->slide_fh(i,j) = 0.0;
}

void sandslide::amr_apply(lexer *p, sediment_fdm *s)
{
    int i,j;
    slice *f = amr_fh();

    if(f==nullptr)
    return;

    SEDSLICELOOP
    {
    s->slide_fh(i,j) += (*f)(i,j);
    s->bedzh(i,j) += (*f)(i,j);
    }
}
