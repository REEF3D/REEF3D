/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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
--------------------------------------------------------------------*/

#include"driver.h"
#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"

void driver::fill_vel(lexer* p, fdm *a, ghostcell *pgc)
{   
    /*
    pgc->dgcpol(p,a->u,p->dgc1,p->dgc1_count,11);
	pgc->dgcpol(p,a->v,p->dgc2,p->dgc2_count,12);
	pgc->dgcpol(p,a->w,p->dgc3,p->dgc3_count,13);
	pgc->dgcpol(p,a->press,p->dgc4,p->dgc4_count,14);
	pgc->dgcpol(p,a->eddyv,p->dgc4,p->dgc4_count,14);
	pgc->dgcpol(p,a->phi,p->dgc4,p->dgc4_count,14);
	pgc->dgcpol(p,a->ro,p->dgc4,p->dgc4_count,14);
	pgc->dgcpol(p,a->visc,p->dgc4,p->dgc4_count,14);
	pgc->dgcpol(p,a->conc,p->dgc4,p->dgc4_count,14);
    pgc->dgcpol(p,a->test,p->dgc4,p->dgc4_count,14);
	
	a->u.ggcpol(p);
	a->v.ggcpol(p);
	a->w.ggcpol(p);
	a->press.ggcpol(p);
	a->eddyv.ggcpol(p);
	a->phi.ggcpol(p);
	a->conc.ggcpol(p);
	a->ro.ggcpol(p);
	a->visc.ggcpol(p);
	a->phi.ggcpol(p);
	a->fb.ggcpol(p);
    a->test.ggcpol(p);
    a->topo.ggcpol(p);
    

    pgc->gcparacox(p,a->phi,50);
	pgc->gcparacox(p,a->phi,50);

	pgc->gcparacox(p,a->topo,150);
	pgc->gcparacox(p,a->topo,150);
    
    pgc->start1(p,a->u,110);
    pgc->start2(p,a->v,111);
	pgc->start3(p,a->w,112);
	
	pgc->dgcpol(p,a->u,p->dgc1,p->dgc1_count,11);
	pgc->dgcpol(p,a->v,p->dgc2,p->dgc2_count,12);
	pgc->dgcpol(p,a->w,p->dgc3,p->dgc3_count,13);
	
	a->u.ggcpol(p);
	a->v.ggcpol(p);
	a->w.ggcpol(p);
    
    
    
    pgc->start1(p,a->u,114);
    pgc->start2(p,a->v,115);
	pgc->start3(p,a->w,116);

	pgc->dgcpol(p,a->u,p->dgc1,p->dgc1_count,11);
	pgc->dgcpol(p,a->v,p->dgc2,p->dgc2_count,12);
	pgc->dgcpol(p,a->w,p->dgc3,p->dgc3_count,13);
	
	a->u.ggcpol(p);
	a->v.ggcpol(p);
	a->w.ggcpol(p);
    */
}







