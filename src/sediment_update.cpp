/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"sediment_f.h"
#include"lexer.h"
#include"fdm.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"sediment_fdm.h"
#include"reinitopo.h"
#include"vrans.h"

void sediment_f::update_cfd(lexer *p, fdm *a,ghostcell *pgc, ioflow *pflow, reinitopo *preto)
{
    topo_zh_update(p,a,pgc,s);
    preto->start(p,a,pgc,a->topo);
    bedchange_update(p, pgc);
    
    volume_calc(p,a,pgc);
    
    pgc->start1(p,a->u,10);
	pgc->start2(p,a->v,11);
	pgc->start3(p,a->w,12);
    
    if(p->mpirank==0)
    cout<<"Topo: update grid..."<<endl;
    
    
    if(p->S10==1 && p->G3==0)
    pgc->topo_update(p,a);
    
    /*if(p->S10==1 && p->G3==1)
    pgc->solid_forcing_topo_update(p,a);*/
    
    if(p->S10==2)
    pvrans->sed_update(p,a,pgc);
    
    pflow->gcio_update(p,a,pgc);
    
    bedlevel(p,a,pgc); 
    
    active_cfd(p,a,pgc);
	
	pgc->start4(p,a->conc,40);
}

void sediment_f::update_sflow(lexer *p, fdm2D *b, ghostcell *pgc, ioflow *pflow)
{
    SLICELOOP4
    b->topobed(i,j) = s->bedzh(i,j);
    
    SLICELOOP4
    b->bed(i,j) = MAX(b->topobed(i,j),b->solidbed(i,j));
    
    pgc->gcsl_start4(p,b->bed,50);
    pgc->gcsl_start4(p,b->topobed,50);
    
    bedchange_update(p, pgc);
    
    active_sflow(p,b,pgc);
    
}