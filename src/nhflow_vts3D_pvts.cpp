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

#include"nhflow_vts3D.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment.h"
#include"nhflow_vel_probe_theory.h"


void nhflow_vts3D::pvts(lexer *p, ghostcell* pgc, sediment *psed)
{	
	int num=0;
    
    if(p->P15==1)
    num = printcount;

    if(p->P15==2)
    num = p->count;
	

	sprintf(name,"./REEF3D_NHFLOW_VTS/REEF3D-NHFLOW-%08i.pvts",num);


	ofstream result;
	result.open(name);

	result<<"<?xml version=\"1.0\"?>"<<endl;
	result<<"<VTKFile type=\"PStructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">"<<endl;
	result<<"<PStructuredGrid WholeExtent=\"0 "<<p->gknox<<" 0 "<<p->gknoy<<" 0 "<<p->gknoz<<"\" GhostLevel=\"0\" Origin=\"0 0 0\" Spacing=\"1 1 1\">"<<endl;
	
    if(p->P16==1)
    {
    result<<"    <FieldData>"<<endl;
    result<<"        <DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\"> "<<p->simtime<<endl;
    result<<"        </DataArray>"<<endl;
    result<<"    </FieldData>"<<endl;
    }
    
	result<<"    <PPointData>"<<endl;
	result<<"        <PDataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\"/>"<<endl;
	result<<"        <PDataArray type=\"Float32\" Name=\"pressure\"/>"<<endl;
    result<<"        <PDataArray type=\"Float32\" Name=\"omega_sig\"/>"<<endl;
    result<<"        <PDataArray type=\"Float32\" Name=\"elevation\"/>"<<endl;
    if(p->P23==1)
	result<<"        <PDataArray type=\"Float32\" Name=\"test\"/>"<<endl;
	result<<"    </PPointData>"<<endl;
	
    result<<"    <PPoints>"<<endl;
	result<<"        <PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>"<<endl;
	result<<"    </PPoints>"<<endl;
    
	for(n=0; n<p->M10; ++n)
	{
    piecename(p,pgc,n);
	extent(p,n);
    result<<"    <Piece Extent=\""<<pextent<<"\" Source=\""<<pname<<"\"/>"<<endl;
	}

	result<<"</PStructuredGrid>"<<endl;
	result<<"</VTKFile>"<<endl;

	result.close();
}

void nhflow_vts3D::piecename(lexer *p, ghostcell *pgc, int n)
{
    int num=0;

    if(p->P15==1)
    num = printcount;

    if(p->P15==2)
    num = p->count;

	sprintf(pname,"REEF3D-NHFLOW-%08i-%06i.vts",num,n+1);
}

void nhflow_vts3D::extent(lexer* p, int n)
{
	sprintf(pextent,"%i %i %i %i %i %i",piextent[0+6*n],piextent[1+6*n],piextent[2+6*n],piextent[3+6*n],piextent[4+6*n],piextent[5+6*n]);
}
void nhflow_vts3D::fextent(lexer* p)
{
	iextent[0]=p->origin_i;
	iextent[1]=p->origin_i+p->knox;
	iextent[2]=p->origin_j;
	iextent[3]=p->origin_j+p->knoy;
	iextent[4]=p->origin_k;
	iextent[5]=p->origin_k+p->knoz;
}