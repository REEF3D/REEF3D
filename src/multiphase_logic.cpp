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

#include"multiphase_f.h"
#include"freesurface_header.h"
#include"multiphase_fluid_update_f.h"
#include"multiphase_fluid_update_rheology.h"
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
	pfsf1 = new levelset_void(p,a,pgc,pheat,pconc);
	pfsf2 = new levelset_void(p,a,pgc,pheat,pconc);
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
	
    
	if(p->F310==0)
	preini = new reini_void(p);
	
	if(p->F310==3)
	preini = new reinifluid_RK3(p,1);

	
	if(p->F310==11 || p->F310==13 || p->F310==14)
	preini = new directreini(p,a);
	

	if(p->F31==0)
	ppls = new particle_pls_void();

	if(p->F31==1 || p->F31==2)
	ppls = new particle_pls(p,a,pgc);
	
	if(p->W90==0)
	pupdate = new multiphase_fluid_update_f(p,a,pgc);
    
    if(p->W90>0)
	pupdate = new multiphase_fluid_update_rheology(p,a,pgc);
}