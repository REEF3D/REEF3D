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

#include"nhflow_vtp_fsf.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"


void nhflow_vtp_fsf::pvtp(lexer *p, fdm_nhf *d, ghostcell* pgc)
{	
	int num=0;

    if(p->P15==1)
    num = printcount;

    if(p->P15==2)
    num = p->count;
	
	sprintf(name,"./REEF3D_NHFLOW_VTP_FSF/REEF3D-NHFLOW-FSF-%08i.pvtp",num);


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
	result<<"<PDataArray type=\"Float64\" NumberOfComponents=\"3\"/>"<<endl;
	result<<"</PPoints>"<<endl;
	
	result<<"<PPointData>"<<endl;
	result<<"<PDataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\"/>"<<endl;
	result<<"<PDataArray type=\"Float32\" Name=\"eta\"/>"<<endl;
    result<<"<PDataArray type=\"Float32\" Name=\"detadt\"/>"<<endl;
	result<<"<PDataArray type=\"Float32\" Name=\"WL\"/>"<<endl;
    result<<"<PDataArray type=\"Float32\" Name=\"breaking\"/>"<<endl;
    result<<"<PDataArray type=\"Float32\" Name=\"coastline\"/>"<<endl;
    result<<"<PDataArray type=\"Float32\" Name=\"wetdry\"/>"<<endl;
    if(p->P23==1)
    result<<"<PDataArray type=\"Float32\" Name=\"test\"/>"<<endl;
    if(p->P28==1)
    result<<"<PDataArray type=\"Float32\" Name=\"fb\"/>"<<endl;
    if(p->P110==1)
    result<<"<PDataArray type=\"Float32\" Name=\"Hs\"/>"<<endl;
    if(p->P131==1)
    result<<"<PDataArray type=\"Float32\" Name=\"wetdry_max\"/>"<<endl;
	result<<"</PPointData>"<<endl;
	
	result<<"<Polys>"<<endl;
    result<<"<DataArray type=\"Int32\"  Name=\"connectivity\"/>"<<endl;
	result<<"<DataArray type=\"Int32\"  Name=\"offsets\" />"<<endl;
    result<<"<DataArray type=\"Int32\"  Name=\"types\" />"<<endl;
	result<<"</Polys>"<<endl;

	for(n=0; n<p->M10; ++n)
	{
    piecename(p,d,pgc,n);
    result<<"<Piece Source=\""<<pname<<"\"/>"<<endl;
	}

	result<<"</PPolyData>"<<endl;
	result<<"</VTKFile>"<<endl;

	result.close();

}

void nhflow_vtp_fsf::piecename(lexer *p, fdm_nhf *d, ghostcell *pgc, int n)
{
    int num=0;


    if(p->P15==1)
    num = printcount;

    if(p->P15==2)
    num = p->count;


    sprintf(pname,"REEF3D-NHFLOW-FSF-%08i-%06i.vtp",num,n+1);

}
