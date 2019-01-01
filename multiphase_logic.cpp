/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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

#include"multiphase_f.h"
#include"freesurface_header.h"
#include"fluid_update_multiphase.h"
#include"lexer.h"
#include"heat.h"
#include"fdm.h"
#include"ghostcell.h"
#include"heat_void.h"
#include"concentration.h"
#include"concentration_void.h"

void multiphase_f::logic(lexer *p, fdm *a, ghostcell *pgc)
{
	pheat = new heat_void(p,a,pgc);
	pconc = new concentration_void(p,a,pgc);
	
	// Free Surface
    if(p->F300==0)
	{
	pfsf1 = new levelset_void(p,a,pgc);
	pfsf2 = new levelset_void(p,a,pgc);
	}

	if(p->F300==1)
	{
	pfsf1 = new levelset_AB2(p,a,pgc,pheat,pconc);
	pfsf2 = new levelset_AB2(p,a,pgc,pheat,pconc);
	}

	if(p->F300==2)
	{
	pfsf1 = new levelset_RK2(p,a,pgc,pheat,pconc);
	pfsf2 = new levelset_RK2(p,a,pgc,pheat,pconc);
	}

	if(p->F300==3)
	{
	pfsf1 = new levelset_RK3(p,a,pgc,pheat,pconc);
	pfsf2 = new levelset_RK3(p,a,pgc,pheat,pconc);
	}
	
	if(p->F300==4)
	{
	pfsf1 = new levelset_RK4(p,a,pgc,pheat,pconc);
	pfsf2 = new levelset_RK4(p,a,pgc,pheat,pconc);
	}

	if(p->F300==11)
	{
	pfsf1 = new levelset_IM1(p,a,pgc,pheat,pconc,ls1);
	pfsf2 = new levelset_IM1(p,a,pgc,pheat,pconc,ls2);
	}

	if(p->F300==12)
	{
	pfsf1 = new levelset_IM2(p,a,pgc,pheat,pconc,ls1);
	pfsf2 = new levelset_IM2(p,a,pgc,pheat,pconc,ls2);
	}
	
	
	if(p->F310==0)
	preini = new reini_void(p);
	
	if(p->F310==1)
	preini = new reini_AB2(p,a);
	
	if(p->F310==3)
	preini = new reini_RK3(p,1);
	
	if(p->F310==4)
	preini = new reini_RK4(p,a);
	
	if(p->F310==5)
	preini = new reinivc_RK3(p);
	
	if(p->F310==7)
	preini = new reinigc_RK3(p,a);
	
	if(p->F310==8)
	preini = new reinigc_RK4(p,a);

	if(p->F310==11 || p->F310==13 || p->F310==14)
	preini = new directreini(p,a);
	
	
	if(p->F31==0)
	ppart = new particle_void();

	if(p->F31==1 || p->F31==2)
	ppart = new particle(p,a,pgc);
	
	
	pupdate = new fluid_update_multiphase(p,a,pgc);
	
	
}