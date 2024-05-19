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

#include"driver.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"mgcslice1.h"
#include"mgcslice2.h"
#include"mgcslice4.h"

void driver::makegrid2D(lexer *p, ghostcell *pgc)
{   
    pgc->gcsl_tpflag(p);
    
    pgc->gcslflagx(p,p->flagslice4);
    
    mgcslice1 m1(p);
	mgcslice2 m2(p);
	mgcslice4 m4(p);
	
	m1.makemgc(p);
    pgc->gcslflagx(p,p->flagslice1);
    m1.gcb_seed(p);
    m1.mgcsetup(p);
    m1.fillmgc(p);
    m1.gcdirfill(p);

    m2.makemgc(p);
    pgc->gcslflagx(p,p->flagslice2);
    m2.gcb_seed(p);
    m2.mgcsetup(p);
    m2.fillmgc(p);
    m2.gcdirfill(p);

    m4.makemgc(p);
    m4.gcb_seed(p);
    m4.mgcsetup(p);
    m4.fillmgc(p);
    m4.gcdirfill(p);
	
	m1.make_ggc(p);
    m1.fill_ggc(p);
	m2.make_ggc(p);
    m2.fill_ggc(p);
    m4.make_ggc(p);
    m4.fill_ggc(p);
    
    pgc->gcsl_setbc1(p);
    pgc->gcsl_setbc2(p);
    pgc->gcsl_setbc4(p);
    
    pgc->gcsl_setbcio(p);
    
    pgc->dgcslini1(p);
	pgc->dgcslini2(p);
	pgc->dgcslini4(p);
}
 
void driver::makegrid2D_cds(lexer *p, ghostcell *pgc, fdm2D *b)
{      
    p->flagini2D();
    p->gridini2D();	

    pgc->sizeS_update(p);
    
    pgc->gcxslupdate(p); 
}
