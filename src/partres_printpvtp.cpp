/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"partres.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<sys/stat.h>
#include<sys/types.h>

void partres::pvtp(lexer* p)
{

    num = printcount;
	
	sprintf(name,"./REEF3D_CFD_SedPart/REEF3D-SedPart-%08i.pvtp",num);
	

	ofstream result;
	result.open(name);

	result<<"<?xml version=\"1.0\"?>\n";
	result<<"<VTKFile type=\"PPolyData\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
	result<<"<PPolyData GhostLevel=\"0\">\n";

	result<<"<FieldData>\n";
	if(p->P16==1)
    {
	result<<"<DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\"> "<<p->simtime;
    result<<"</DataArray>\n";
	}
	result<<"</FieldData>\n";

	result<<"<PPointData>\n";
	result<<"<PDataArray type=\"Float32\" Name=\"Flag\"/>\n";
    if(p->P23==1)
    result<<"<PDataArray type=\"Float32\" Name=\"Test\"/>\n";
	result<<"<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\"/>\n";
	result<<"<PDataArray type=\"Float32\" Name=\"radius\"/>\n";
    result<<"<DataArray type=\"Float32\" Name=\"fluid velocity\" NumberOfComponents=\"3\"/>\n";
	result<<"<PDataArray type=\"Float32\" Name=\"bedChange\"/>\n";
	result<<"</PPointData>\n";

	result<<"<PPoints>\n";
	result<<"<PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>\n";
	result<<"</PPoints>\n";

	for(int n=0; n<p->M10; ++n)
	{
    piecename_pos(p,n);
    result<<"<Piece Source=\""<<pname<<"\"/>\n";
	}

	result<<"</PPolyData>\n";
	result<<"</VTKFile>"<<std::flush;

	result.close();
}

void partres::piecename_pos(lexer* p, int n)
{
	sprintf(pname,"REEF3D-SedPart-%08i-%06i.vtp",printcount,n+1);
}

/// @brief Setting name of indivdual vtp file
void partres::header_pos(lexer* p)
{
	sprintf(name,"./REEF3D_CFD_SedPart/REEF3D-SedPart-%08i-%06i.vtp",printcount,p->mpirank+1);	
}
