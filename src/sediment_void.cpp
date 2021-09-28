/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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

#include"sediment_void.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
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

void sediment_void::start(lexer *p, fdm *a, convection *pconvec, ghostcell *pgc, ioflow *pflow,
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

void sediment_void::print_3D_bedshear(lexer* p, fdm *a, ghostcell *pgc, ofstream &result)
{	
}

void sediment_void::name_pvtu_bedshear(lexer *p, fdm *a, ghostcell *pgc, ofstream &result)
{
}

void sediment_void::name_vtu_bedshear(lexer *p, fdm *a, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
}

void sediment_void::offset_vtu_bedshear(lexer *p, fdm *a, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
}

void sediment_void::print_3D_parameter1(lexer* p, fdm *a, ghostcell *pgc, ofstream &result)
{	
}

void sediment_void::name_pvtu_parameter1(lexer *p, fdm *a, ghostcell *pgc, ofstream &result)
{
}

void sediment_void::name_vtu_parameter1(lexer *p, fdm *a, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
}

void sediment_void::offset_vtu_parameter1(lexer *p, fdm *a, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
}

void sediment_void::print_3D_parameter2(lexer* p, fdm *a, ghostcell *pgc, ofstream &result)
{	
}

void sediment_void::name_pvtu_parameter2(lexer *p, fdm *a, ghostcell *pgc, ofstream &result)
{
}

void sediment_void::name_vtu_parameter2(lexer *p, fdm *a, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
}

void sediment_void::offset_vtu_parameter2(lexer *p, fdm *a, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
}

double sediment_void::bedshear_point(lexer *p, fdm *a,ghostcell *pgc)
{
	return 0.0;
}
