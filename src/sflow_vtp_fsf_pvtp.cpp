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

#include"sflow_vtp_fsf.h"
#include"sflow_turbulence.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"sediment.h"

void sflow_vtp_fsf::pvtp(lexer *p, fdm2D* b, ghostcell* pgc, sflow_turbulence *pturb, sediment *psed)
{	
	int num=0;

    if(p->P15==1)
    num = p->printcount;

    if(p->P15==2)
    num = p->count;
	
	sprintf(name,"./REEF3D_SFLOW_VTP_FSF/REEF3D-SFLOW-FSF-%08i.pvtp",num);

	ofstream result;
	result.open(name);

	result<<"<?xml version=\"1.0\"?>"<<endl;
	result<<"<VTKFile type=\"PPolyData\" version=\"0.1\" byte_order=\"LittleEndian\">"<<endl;
	result<<"<PPolyData  GhostLevel=\"0\">"<<endl;
	
	result<<"<PPoints>"<<endl;
	result<<"<PDataArray type=\"Float64\" NumberOfComponents=\"3\"/>"<<endl;
	result<<"</PPoints>"<<endl;
	
	result<<"<PPointData>"<<endl;
	result<<"<PDataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\"/>"<<endl;
	result<<"<PDataArray type=\"Float32\" Name=\"pressure\"/>"<<endl;
	result<<"<PDataArray type=\"Float32\" Name=\"eddyv\"/>"<<endl;
    pturb->name_pvtp(p,b,pgc,result);
	result<<"<PDataArray type=\"Float32\" Name=\"eta\"/>"<<endl;
	result<<"<PDataArray type=\"Float32\" Name=\"waterlevel\"/>"<<endl;
    result<<"<PDataArray type=\"Float32\" Name=\"wetdry\"/>"<<endl;
    result<<"<PDataArray type=\"Float32\" Name=\"breaking\"/>"<<endl;
    
    if(p->P23==1)
    result<<"<PDataArray type=\"Float32\" Name=\"test\"/>"<<endl;
    if(p->P28==1)
    result<<"<PDataArray type=\"Float32\" Name=\"fb\"/>"<<endl;
    if(p->P110==1)
    result<<"<PDataArray type=\"Float32\" Name=\"Hs\"/>"<<endl;
    
	result<<"</PPointData>"<<endl;
	
	result<<"<Polys>"<<endl;
    result<<"<DataArray type=\"Int32\"  Name=\"connectivity\"/>"<<endl;
	result<<"<DataArray type=\"Int32\"  Name=\"offsets\" />"<<endl;
    result<<"<DataArray type=\"Int32\"  Name=\"types\" />"<<endl;
	result<<"</Polys>"<<endl;

	for(n=0; n<p->M10; ++n)
	{
    piecename(p,b,pgc,n);
    result<<"<Piece Source=\""<<pname<<"\"/>"<<endl;
	}

	result<<"</PPolyData>"<<endl;
	result<<"</VTKFile>"<<endl;

	result.close();

}

void sflow_vtp_fsf::piecename(lexer *p, fdm2D *b, ghostcell *pgc, int n)
{
    int num=0;


    if(p->P15==1)
    num = p->printcount;

    if(p->P15==2)
    num = p->count;

	sprintf(pname,"REEF3D-SFLOW-FSF-%08i-%06i.vtp",num,n+1);
}
