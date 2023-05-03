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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"sflow_turb_io_void.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"

sflow_turb_io_void::sflow_turb_io_void(lexer* p)
{

}

sflow_turb_io_void::~sflow_turb_io_void()
{
}

void sflow_turb_io_void::print_2D(lexer *p, fdm2D *b, ghostcell *pgc, ofstream &result)
{

}

void sflow_turb_io_void::kinget(int ii, int jj, double val)
{

}

void sflow_turb_io_void::epsget(int ii, int jj, double val)
{
    
}
    
double sflow_turb_io_void::kinval(int ii, int jj)
{

}
    
double sflow_turb_io_void::epsval(int ii, int jj)
{

}
    
void sflow_turb_io_void::name_pvtp(lexer *p, fdm2D *b, ghostcell *pgc,ofstream &result)
{

}

void sflow_turb_io_void::name_vtp(lexer *p, fdm2D *b, ghostcell *pgc,ofstream &result, int *offset, int &n)
{

}
    
void sflow_turb_io_void::offset_vtp(lexer *p, fdm2D *b, ghostcell *pgc,ofstream &result, int *offset, int &n)
{

}