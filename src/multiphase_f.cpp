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
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"solver.h"
#include"ioflow.h"
#include"reini.h"
#include"picard.h"
#include"multiphase_fluid_update.h"
#include"picard_f.h"
#include"picard_lsm.h"
#include"picard_void.h"
#include"print_wsf.h"

multiphase_f::multiphase_f(lexer* p, fdm *a, ghostcell* pgc) : ls1(p), ls2(p)
{
	logic(p,a,pgc);
	
	pwsf1=new print_wsf(p,a,pgc,1);
	pwsf2=new print_wsf(p,a,pgc,2);
}

multiphase_f::~multiphase_f()
{
}

void multiphase_f::start(lexer *p, fdm *a, ghostcell *pgc, convection *pmpconvec, solver *psolv, ioflow *pflow, reini* preini2, particle_corr* ppls, printer *pprint)
{
	pfsf1->start(a,p,pmpconvec,psolv,pgc,pflow,preini,ppls,ls1);
	pfsf2->start(a,p,pmpconvec,psolv,pgc,pflow,preini,ppls,ls2);	
	
	update(p,a,pgc);
}

