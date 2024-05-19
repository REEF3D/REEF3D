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

#include"print_averaging_v.h"
#include"lexer.h"
#include"fdm.h"
#include"heat.h"
#include"ghostcell.h"
#include<iostream>
#include<fstream>
#include<sys/stat.h>
#include<sys/types.h>

print_averaging_v::print_averaging_v(lexer *p, fdm* a, ghostcell *pgc) 
{

}

print_averaging_v::~print_averaging_v()
{

}

void print_averaging_v::averaging(lexer *p, fdm *a, ghostcell *pgc, heat *pheat)
{

}

void print_averaging_v::offset_vtu(lexer *p, fdm *a, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
        
}

void print_averaging_v::name_vtu(lexer *p, fdm *a, ghostcell *pgc, ofstream &result, int *offset, int &n)
{

}

void print_averaging_v::name_pvtu(lexer *p, fdm *a, ghostcell *pgc, ofstream &result)
{

}

void print_averaging_v::print_3D(lexer* p, fdm *a, ghostcell *pgc, ofstream &result)
{

}

