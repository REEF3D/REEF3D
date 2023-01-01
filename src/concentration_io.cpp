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

#include"concentration_io.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

concentration_io::concentration_io(lexer *p, fdm *a) : C(p)
{
}

concentration_io::~concentration_io()
{
}

void concentration_io::print_3D(lexer* p, fdm *a, ghostcell *pgc, ofstream &result)
{
	
	iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));
	
	TPLOOP
	{
	ffn=float(p->ipol4(C));
	result.write((char*)&ffn, sizeof (float));
	}	
	
	iin=4*(p->pointnum);
    result.write((char*)&iin, sizeof (int));

	TPLOOP
	{
	ffn=float(p->ipol4(a->ro));
	result.write((char*)&ffn, sizeof (float));
	}
}

double concentration_io::val(int ii, int jj, int kk)
{
    double val;

    val=C(ii,jj,kk);

    return val;
}

void concentration_io::name_pvtu(lexer *p, fdm *a, ghostcell *pgc, ofstream &result)
{
    result<<"<PDataArray type=\"Float32\" Name=\"C\"/>"<<endl;
	result<<"<PDataArray type=\"Float32\" Name=\"rho\"/>"<<endl;
}

void concentration_io::name_vtu(lexer *p, fdm *a, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
    result<<"<DataArray type=\"Float32\" Name=\"C\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	result<<"<DataArray type=\"Float32\" Name=\"rho\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
}

void concentration_io::offset_vtu(lexer *p, fdm *a, ghostcell *pgc, ofstream &result, int *offset, int &n)
{
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
	offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
}
