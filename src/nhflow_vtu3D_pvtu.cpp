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

#include"nhflow_vtu3D.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment.h"
#include"nhflow_turbulence.h"

void nhflow_vtu3D::pvtu(lexer *p, fdm_nhf *d, ghostcell* pgc, nhflow_turbulence *pnhfturb, sediment *psed)
{	
	int num=0;
    
    if(p->P15==1)
    num = printcount;

    if(p->P15==2)
    num = p->count;
	

	sprintf(name,"./REEF3D_NHFLOW_VTU/REEF3D-NHFLOW-%08i.pvtu",num);


	ofstream result;
	result.open(name);

	result<<"<?xml version=\"1.0\"?>"<<endl;
	result<<"<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"<<endl;
	result<<"<PUnstructuredGrid GhostLevel=\"0\">"<<endl;
	
    if(p->P16==1)
    {
    result<<"<FieldData>"<<endl;
    result<<"<DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\"> "<<p->simtime<<endl;
    result<<"</DataArray>"<<endl;
    result<<"</FieldData>"<<endl;
    }
    
	result<<"<PPointData>"<<endl;
	result<<"<PDataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\"/>"<<endl;
	result<<"<PDataArray type=\"Float32\" Name=\"pressure\"/>"<<endl;
    pnhfturb->name_pvtu(p,d,pgc,result);
    result<<"<PDataArray type=\"Float32\" Name=\"omega_sig\"/>"<<endl;
    result<<"<PDataArray type=\"Float32\" Name=\"elevation\"/>"<<endl;
    if(p->P23==1)
	result<<"<PDataArray type=\"Float32\" Name=\"test\"/>"<<endl;
    if(p->P110==1)
	result<<"<PDataArray type=\"Float32\" Name=\"Hs\"/>"<<endl;
    if(p->P25==1)
	result<<"<PDataArray type=\"Float32\" Name=\"solid\"/>"<<endl;
    if(p->P25==1 || p->P28==1)
    result<<"<PDataArray type=\"Float32\" Name=\"Heaviside\"/>"<<endl;
    if(p->P28==1)
	result<<"<PDataArray type=\"Float32\" Name=\"floating\"/>"<<endl;
	result<<"</PPointData>"<<endl;
	
    result<<"<PPoints>"<<endl;
	result<<"<PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>"<<endl;
	result<<"</PPoints>"<<endl;
    
	result<<"<Cells>"<<endl;
    result<<"<DataArray type=\"Int32\"  Name=\"connectivity\"/>"<<endl;
	result<<"<DataArray type=\"Int32\"  Name=\"offsets\" />"<<endl;
    result<<"<DataArray type=\"Int32\"  Name=\"types\" />"<<endl;
	result<<"</Cells>"<<endl;

	for(n=0; n<p->M10; ++n)
	{
    piecename(p,pgc,n);
    result<<"<Piece Source=\""<<pname<<"\"/>"<<endl;
	}

	result<<"</PUnstructuredGrid>"<<endl;
	result<<"</VTKFile>"<<endl;

	result.close();
}

void nhflow_vtu3D::piecename(lexer *p, ghostcell *pgc, int n)
{
    int num=0;

    if(p->P15==1)
    num = printcount;

    if(p->P15==2)
    num = p->count;

	sprintf(pname,"REEF3D-NHFLOW-%08i-%06i.vtu",num,n+1);
}
