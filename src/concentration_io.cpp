/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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
#include<cstring>

concentration_io::concentration_io(lexer *p, fdm *a) : C(p)
{
}

concentration_io::~concentration_io()
{
}

void concentration_io::print_3D(lexer* p, fdm *a, ghostcell *pgc, std::vector<char> &buffer, size_t &m)
{
	
	iin=4*(p->pointnum);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
	
	TPLOOP
	{
	ffn=float(p->ipol4(C));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
	m+=sizeof(float);
	}	
	
	iin=4*(p->pointnum);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);

	TPLOOP
	{
	ffn=float(p->ipol4(a->ro));
	std::memcpy(&buffer[m],&ffn,sizeof(float));
	m+=sizeof(float);
	}
}

double concentration_io::val(int ii, int jj, int kk)
{
    double val;

    val=C(ii,jj,kk);

    return val;
}

void concentration_io::name_ParaView_parallel(lexer *p, ofstream &result)
{
    result<<"<PDataArray type=\"Float32\" Name=\"C\"/>\n";
	result<<"<PDataArray type=\"Float32\" Name=\"rho\"/>\n";
}

void concentration_io::name_ParaView(lexer *p, ostream &result, int *offset, int &n)
{
    result<<"<DataArray type=\"Float32\" Name=\"C\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
	result<<"<DataArray type=\"Float32\" Name=\"rho\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
}

void concentration_io::offset_ParaView(lexer *p, int *offset, int &n)
{
    offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
	offset[n]=offset[n-1]+4*(p->pointnum)+4;
	++n;
}
