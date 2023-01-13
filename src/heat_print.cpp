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

#include"heat_print.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

heat_print::heat_print(lexer *p, fdm *a) : T(p)
{
}

heat_print::~heat_print()
{
}

void heat_print::print_3D(lexer* p, fdm *a, ghostcell *pgc, ofstream &result)
{
	
	iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));
	
	TPLOOP
	{
	ffn=float(p->ipol4_a(T));
	result.write((char*)&ffn, sizeof (float));
	}
}

double heat_print::val(int ii, int jj, int kk)
{
    double val;

    val=T(ii,jj,kk);

    return val;
}

void heat_print::name_pvtu(lexer *p, fdm *a, ghostcell *pgc, ofstream &result)
{
    result<<"<PDataArray type=\"Float32\" Name=\"T\"/>"<<endl;
}

void heat_print::name_vtu(lexer *p, fdm *a, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
    result<<"<DataArray type=\"Float32\" Name=\"T\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
}

void heat_print::offset_vtu(lexer *p, fdm *a, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
}
