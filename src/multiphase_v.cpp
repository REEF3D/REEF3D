/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"multiphase_v.h"
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


multiphase_v::multiphase_v()
{
}

multiphase_v::~multiphase_v()
{
}

void multiphase_v::ini(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow, printer *pprint, convection *pconvec, solver *psolv)
{
}

void multiphase_v::start(lexer *p, fdm *a, ghostcell *pgc, convection *pconvec, solver *psolv, ioflow *pflow, reini* preini, particle_corr* ppls, printer *pprint)
{
}

void multiphase_v::update(lexer *p, fdm *a, ghostcell *pgc)
{
}

void multiphase_v::print_3D(lexer*, fdm*, ghostcell*,ofstream&)
{	
}

void multiphase_v::print_file(lexer *p, fdm *a, ghostcell *pgc)
{
}

void multiphase_v::nodefill(lexer*,fdm*,ghostcell*,field&)
{
	
}

double multiphase_v::ls1val(int,int,int)
{
	double val=0.0;
	
	return val;
}

double multiphase_v::ls2val(int,int,int)
{
	double val=0.0;
	
	return val;
}

double multiphase_v::ccipol_ls1val(lexer*,ghostcell*,double,double,double)
{
	double val=0.0;
	
	return val;
}

double multiphase_v::ccipol_ls2val(lexer*,ghostcell*,double,double,double)
{
	double val=0.0;
	
	return val;
}

void multiphase_v::ls1get(int,int,int,double)
{
	
}

void multiphase_v::ls2get(int,int,int,double)
{
	
}


void multiphase_v::name_pvtu(lexer*, fdm*, ghostcell*,ofstream&)
{
	
}

void multiphase_v::name_vtu(lexer*, fdm*, ghostcell*,ofstream&, int*, int &)
{
	
}

void multiphase_v::offset_vtu(lexer*, fdm*, ghostcell*,ofstream&, int*, int &)
{
	
}

