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

#include"fsf_vtp.h"
#include<string>
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void fsf_vtp::pvtp(lexer* p, fdm* a, ghostcell* pgc)
{
    int num=0;

    if(p->P15==1)
    num = fsfprintcount;

    if(p->P15==2)
    num = p->count;
	
	sprintf(name,"./REEF3D_CFD_FSF/REEF3D-CFD-FSF-%08i.pvtp",num);

	ofstream result;
	result.open(name);

	result<<"<?xml version=\"1.0\"?>"<<endl;
	result<<"<VTKFile type=\"PPolyData\" version=\"0.1\" byte_order=\"LittleEndian\">"<<endl;
	result<<"<PPolyData  GhostLevel=\"0\">"<<endl;
    
    if(p->P16==1)
    {
    result<<"<FieldData>"<<endl;
    result<<"<DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\"> "<<p->simtime<<endl;
    result<<"</DataArray>"<<endl;
    result<<"</FieldData>"<<endl;
    }

	result<<"<PPoints>"<<endl;
	result<<"<PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>"<<endl;
	result<<"</PPoints>"<<endl;
	
	result<<"<PPointData>"<<endl;
	result<<"<PDataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\"/>"<<endl;
	result<<"<PDataArray type=\"Float32\" Name=\"elevation\"/>"<<endl;
	result<<"</PPointData>"<<endl;
	
	result<<"<Polys>"<<endl;
    result<<"<DataArray type=\"Int32\"  Name=\"connectivity\"/>"<<endl;
    ++n;
	result<<"<DataArray type=\"Int32\"  Name=\"offsets\"/>"<<endl;
	++n;
    result<<"<DataArray type=\"Int32\"  Name=\"types\"/>"<<endl;
	result<<"</Polys>"<<endl;

	for(n=0; n<p->M10; ++n)
	{
    piecename(p,a,pgc,n);
    result<<"<Piece Source=\""<<pname<<"\"/>"<<endl;
	}

	result<<"</PPolyData>"<<endl;
	result<<"</VTKFile>"<<endl;

	result.close();
}

void fsf_vtp::piecename(lexer* p, fdm* a,  ghostcell* pgc, int n)
{
    int num=0;


    if(p->P15==1)
    num = fsfprintcount;

    if(p->P15==2)
    num = p->count;

	sprintf(pname,"REEF3D-CFD-FSF-%08i-%06i.vtp",num,n+1);

}
