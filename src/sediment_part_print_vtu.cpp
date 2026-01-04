/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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
Author: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#include"sediment_part.h"
#include"lexer.h"
#include"CPM.h"
#include"sediment_fdm.h"

void sediment_part::name_ParaView_parallel_CPM(lexer *p, ofstream &result)
{
    pst->name_ParaView_parallel_CPM(p,result);
}

void sediment_part::name_ParaView_CPM(lexer *p, ostream &result, int *offset, int &n)
{
    pst->name_ParaView_CPM(p,result,offset,n);
}

void sediment_part::offset_ParaView_CPM(lexer *p, int *offset, int &n)
{
    pst->offset_ParaView_CPM(p,offset,n);
}

void sediment_part::print_3D_CPM(lexer* p, ghostcell *pgc, vector<char> &buffer, size_t &m)
{	
    pst->print_3D_CPM(p,pgc,buffer,m);
}