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
#include"fdm.h"
#include"ghostcell.h"
#include"mgc1.h"
#include"mgc2.h"
#include"mgc3.h"
#include"mgc4.h"
#include"mgc4a.h"
#include"mgc6.h"
#include"cart1.h"
#include"cart2.h"
#include"cart3.h"
#include"cart4.h"
#include"cart4a.h"

void driver::makegrid(lexer *p, ghostcell *pgc)
{	
	mgc1 m1(p);
	mgc2 m2(p);
	mgc3 m3(p);
	mgc4 m4(p);
	mgc4a m4a(p);
    mgc6 m6(p);
    
	pgc->flagx(p,p->flag1);
    pgc->flagx(p,p->flag2);
    pgc->flagx(p,p->flag3);
    pgc->flagx(p,p->flag4);
    pgc->flagx(p,p->flag);
	pgc->gcxupdate(p);
    

	m1.makemgc(p);
    m1.fillgcb(p);
    m1.extragcb(p);
    m1.mgcsetup(p);
    m1.fillmgc(p);
    m1.gcdirfill(p);
    
    m2.makemgc(p);
    m2.fillgcb(p);
    m2.extragcb(p);
    m2.mgcsetup(p);
    m2.fillmgc(p);
    m2.gcdirfill(p);
    
    m3.makemgc(p);
    m3.fillgcb(p);
    m3.extragcb(p);
    m3.mgcsetup(p);
    m3.fillmgc(p);
    m3.gcdirfill(p);
    
    m4.makemgc(p);
    m4.mgcsetup(p);
    m4.fillmgc(p);
    m4.gcdirfill(p);
	m4.gcsidefill(p);
    
    m4a.makemgc(p);
    m4a.fillgcb(p);
    m4a.mgcsetup(p);
    m4a.fillmgc(p);
    m4a.gcdirfill(p);
    
    m6.makemgc(p);
    m6.mgcsetup(p);
    m6.fillmgc(p);
    m6.gcdirfill(p);
    
	m1.make_ggc(p);
    m1.fill_ggc(p);
	m2.make_ggc(p);
    m2.fill_ggc(p);
	m3.make_ggc(p);
    m3.fill_ggc(p);
    m4.make_ggc(p);
    m4.fill_ggc(p);
    m4a.make_ggc(p);
    m4a.fill_ggc(p);
    
    m1.make_dgc(p);
    m2.make_dgc(p);
    m3.make_dgc(p);
    m4.make_dgc(p);
    
    m1.fill_dgc(p);
    m2.fill_dgc(p);
    m3.fill_dgc(p);
    m4.fill_dgc(p);
    
    p->vecsize(pgc);
}
	
void driver::makegrid_cds()
{	
	pgc->sizeM_update(p,a);
    
    pgc->column_pt4_update(p,a);
    pgc->column_pt4a_update(p,a);
    pgc->column_pt6_update(p,a);
}
	
