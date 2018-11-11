/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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
--------------------------------------------------------------------*/

#include"sediment_void.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"discrete.h"
#include"ioflow.h"
#include"topo.h"
#include"reinitopo.h"
#include"suspended.h"
#include"bedload.h"

sediment_void::sediment_void()
{

}

sediment_void::~sediment_void()
{

}

void sediment_void::start(lexer *p, fdm *a, discrete *pdisc, ghostcell *pgc, ioflow *pflow,
                                        topo *ptopo, reinitopo *preto, suspended *psusp, bedload *pbed)
{

}

void sediment_void::update(lexer *p, fdm *a,ghostcell *pgc, ioflow *pflow)
{
}

void sediment_void::relax(lexer *p, fdm *a,ghostcell *pgc)
{
}

void sediment_void::ini(lexer *p, fdm *a,ghostcell *pgc)
{
}

double sediment_void::bedshear_point(lexer *p, fdm *a, ghostcell *pgc)
{
	double val=0.0;

	return val;	
}


void sediment_void::print_3D(lexer* p, fdm *a, ghostcell *pgc, ofstream &result)
{	
}

void sediment_void::name_pvtu(lexer *p, fdm *a, ghostcell *pgc, ofstream &result)
{
}

void sediment_void::name_vtu(lexer *p, fdm *a, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
}

void sediment_void::offset_vtu(lexer *p, fdm *a, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
}
