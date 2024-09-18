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
for more details->

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"partres2.h"
#include"lexer.h"
#include"sediment_fdm.h"
#include"ghostcell.h"

void partres2::print_particles(lexer* p, sediment_fdm *s)
{
    if((p->count%p->Q181==0 || p->count==0) && (p->Q180==1 && p->Q181>0 && p->Q182<0.0))
	{
    print_vtp(p,s);
	++printcount;
	}
    
    if((p->simtime>p->partprinttime || p->count==0) && (p->Q180==1 && p->Q181<0 && p->Q182>0.0))
    {
    print_vtp(p,s);
    p->partprinttime+=p->Q182;
	++printcount;
    }

	if(((p->count%p->Q181==0 && p->count>=p->Q182) || p->count==0) && p->Q180==2)
	{
	print_vtp(p,s);
	++printcount;
	}
    
}

void partres2::print_vtp(lexer* p, sediment_fdm *s)
{
	int numpt=0;

	for(n=0;n<P.index;++n)
    if(P.Flag[n]>0)
	numpt++;

	//cout<<"PSed-"<<p->mpirank<<"| printed: "<<numpt<<" not printed: "<<P.size-numpt<<" | capcacity: "<<P.capacity<<endl;

	int count;
	int n=0;
    int offset[100];
	int iin;
	float ffn;
	
	if(p->mpirank==0)
	pvtp(p);

    header_pos(p);

	ofstream result;
	result.open(name, ios::binary);


	offset[n]=0;
	++n;
	
	offset[n]=offset[n-1]+4*(numpt)+4; //flag
	++n;
	offset[n]=offset[n-1]+4*(numpt)*3+4; //velocity
	++n;
	offset[n]=offset[n-1]+4*(numpt)+4; //radius
	++n;
    offset[n]=offset[n-1]+4*(numpt)*3+4; //fluid velocity
	++n;
	offset[n]=offset[n-1]+4*(numpt)+4; //bedChange
	++n;
	

    offset[n]=offset[n-1]+4*(numpt)*3+4; //xyz
    ++n;
    offset[n]=offset[n-1]+4*(numpt)+4; //connectivitey
    ++n;
	offset[n]=offset[n-1]+4*(numpt)+4; //offset connectivity
    ++n;

	//---------------------------------------------
	n=0;
	result<<"<?xml version=\"1.0\"?>\n";
	result<<"<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
	result<<"<PolyData>\n";
	result<<"<Piece NumberOfPoints=\""<<numpt<<"\" NumberOfVerts=\""<<numpt<<"\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n";
	
	result<<"<FieldData>\n";
	if(p->P16==1)
    {
	result<<"<DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\"> "<<p->simtime<<endl;
    result<<"</DataArray>\n";
	}
	result<<"</FieldData>\n";
	
	result<<"<PointData >\n";
	result<<"<DataArray type=\"Float32\" Name=\"Flag\"  format=\"appended\" offset=\""<<offset[n]<<"\" />\n";
    ++n;
	result<<"<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset[n]<<"\" />\n";
	++n;
    result<<"<DataArray type=\"Float32\" Name=\"radius\" format=\"appended\" offset=\""<<offset[n]<<"\" />\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"fluid velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset[n]<<"\" />\n";
    ++n;
	result<<"<DataArray type=\"Float32\" Name=\"bedChange\" format=\"appended\" offset=\""<<offset[n]<<"\" />\n";
    ++n;
	result<<"</PointData>\n";
	

    result<<"<Points>\n";
    result<<"<DataArray type=\"Float32\"  NumberOfComponents=\"3\"  format=\"appended\" offset=\""<<offset[n]<<"\" />\n";
    ++n;
    result<<"</Points>\n";
	

    result<<"<Verts>\n";
	result<<"<DataArray type=\"Int32\"  Name=\"connectivity\"  format=\"appended\" offset=\""<<offset[n]<<"\" />\n";
    ++n;
	result<<"<DataArray type=\"Int32\"  Name=\"offsets\"  format=\"appended\" offset=\""<<offset[n]<<"\" />\n";
	++n;
	result<<"</Verts>\n";

    result<<"</Piece>\n";
    result<<"</PolyData>\n";

	//----------------------------------------------------------------------------
    result<<"<AppendedData encoding=\"raw\">\n"<<"_";
	
	// flag
    iin=4*(numpt);
    result.write((char*)&iin, sizeof (int));
	for(n=0;n<P.index;++n)
	if(P.Flag[n]>=0)
	{
		ffn=float(P.Flag[n]);
		result.write((char*)&ffn, sizeof (float));
	}

	// velocities
	iin=4*(numpt)*3;
	result.write((char*)&iin, sizeof (int));
    for(n=0;n<P.index;++n)
    if(P.Flag[n]>=0)
	{
	ffn=float(P.U[n]);
	result.write((char*)&ffn, sizeof (float));

	ffn=float(P.V[n]);
	result.write((char*)&ffn, sizeof (float));

	ffn=float(P.W[n]);
	result.write((char*)&ffn, sizeof (float));
	}

	// radius
    iin=4*(numpt);
    result.write((char*)&iin, sizeof (int));
	for(n=0;n<P.index;++n)
	if(P.Flag[n]>=0)
	{
		ffn=float(P.d50/2);
		result.write((char*)&ffn, sizeof (float));
	}

    // fluid velocities
	iin=4*(numpt)*3;
	result.write((char*)&iin, sizeof (int));
    for(n=0;n<P.index;++n)
    if(P.Flag[n]>=0)
	{
	ffn=float(P.Uf[n]);
	result.write((char*)&ffn, sizeof (float));

	ffn=float(P.Vf[n]);
	result.write((char*)&ffn, sizeof (float));

	ffn=float(P.Wf[n]);
	result.write((char*)&ffn, sizeof (float));
	}

	// bedChange
    iin=4*(numpt);
    result.write((char*)&iin, sizeof (int));
	for(n=0;n<P.index;++n)
	if(P.Flag[n]>=0)
	{
		//ffn=float(p->ccslipol4(s->bedch,P.X[n],P.Y[n]));
         ffn=float(p->ccslipol4(s->bedzh,P.X[n],P.Y[n])-p->ccslipol4(s->bedzh0,P.X[n],P.Y[n]));
		result.write((char*)&ffn, sizeof (float));
	}

	//  XYZ
	iin=4*(numpt)*3;
	result.write((char*)&iin, sizeof (int));
    for(n=0;n<P.index;++n)
    if(P.Flag[n]>=0)
	{
	ffn=float(P.X[n]);
	result.write((char*)&ffn, sizeof (float));

	ffn=float(P.Y[n]);
	result.write((char*)&ffn, sizeof (float));

	ffn=float(P.Z[n]);
	result.write((char*)&ffn, sizeof (float));
	}
	
	//  Connectivity
	count=0;
    iin=4*(numpt);
    result.write((char*)&iin, sizeof (int));
	for(n=0;n<P.index;++n)
	if(P.Flag[n]>=0)
	{
	iin=int(count);
	result.write((char*)&iin, sizeof (int));
	++count;
	}

	//  Offset of Connectivity
	count=1;
    iin=4*(numpt);
    result.write((char*)&iin, sizeof (int));
	for(n=0;n<P.index;++n)
    if(P.Flag[n]>=0)
	{
	iin=int(count);
	result.write((char*)&iin, sizeof (int));
	++count;
	}

	result<<"\n</AppendedData>\n";
    result<<"</VTKFile>"<<flush;

	result.close();
	}
