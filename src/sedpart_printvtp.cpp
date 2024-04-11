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

#include"sedpart.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

/// @brief Printing contol function
void sedpart::print_particles(lexer* p)
{
    if((p->count%p->Q181==0 || p->count==0) && (p->Q180==1 && p->Q181>0 && p->Q182<0.0))
	{
    print_vtp(p);
	++printcount;
	}
    
    if((p->simtime>p->partprinttime || p->count==0) && (p->Q180==1 && p->Q181<0 && p->Q182>0.0))
    {
    print_vtp(p);
    p->partprinttime+=p->Q182;
	++printcount;
    }
    
}

/// @brief Printing particle as vtp
void sedpart::print_vtp(lexer* p)
{
	int numpt=0;
	const int print_flag=p->Q183;

	PARTLOOP
	if(PP.Flag[n]>=print_flag)
	numpt++;

	cout<<"PSed-"<<p->mpirank<<"| printed: "<<numpt<<" not printed: "<<PP.size-numpt<<" | capcaity: "<<PP.capacity<<endl;

	int count;
	int n=0;
    int offset[100];
	int iin;
	float ffn;
	
	if(p->mpirank==0)
	pvtp_pos(p);

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
	

    offset[n]=offset[n-1]+4*(numpt)*3+4; //xyz
    ++n;
    offset[n]=offset[n-1]+4*(numpt)+4; //connectivitey
    ++n;
	offset[n]=offset[n-1]+4*(numpt)+4; //offset connectivity
    ++n;

	//---------------------------------------------
	n=0;
	result<<"<?xml version=\"1.0\"?>"<<endl;
	result<<"<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\">"<<endl;
	result<<"<PolyData>"<<endl;
	result<<"<Piece NumberOfPoints=\""<<numpt<<"\" NumberOfVerts=\""<<numpt<<"\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">"<<endl;
	
	result<<"<FieldData>"<<endl;
	result<<"<DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\"> "<<p->simtime<<endl;
    result<<"</DataArray>"<<endl;
	result<<"</FieldData>"<<endl;
	
	result<<"<PointData >"<<endl;
	result<<"<DataArray type=\"Float32\" Name=\"Flag\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	result<<"<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
	++n;
    result<<"<DataArray type=\"Float32\" Name=\"radius\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	result<<"</PointData>"<<endl;
	

    result<<"<Points>"<<endl;
    result<<"<DataArray type=\"Float32\"  NumberOfComponents=\"3\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"</Points>"<<endl;
	

    result<<"<Verts>"<<endl;
	result<<"<DataArray type=\"Int32\"  Name=\"connectivity\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	result<<"<DataArray type=\"Int32\"  Name=\"offsets\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
	++n;
	result<<"</Verts>"<<endl;

    result<<"</Piece>"<<endl;
    result<<"</PolyData>"<<endl;

	//----------------------------------------------------------------------------
    result<<"<AppendedData encoding=\"raw\">"<<endl<<"_";
	
	// flag
    iin=4*(numpt);
    result.write((char*)&iin, sizeof (int));
	PARTLOOP
	if(PP.Flag[n]>=print_flag)
	{
		ffn=float(PP.Flag[n]);
		result.write((char*)&ffn, sizeof (float));
	}

	// velocities
	iin=4*(numpt)*3;
	result.write((char*)&iin, sizeof (int));
    PARTLOOP
    if(PP.Flag[n]>=print_flag)
	{
	ffn=float(PP.U[n]);
	result.write((char*)&ffn, sizeof (float));

	ffn=float(PP.V[n]);
	result.write((char*)&ffn, sizeof (float));

	ffn=float(PP.W[n]);
	result.write((char*)&ffn, sizeof (float));
	}

	// radius
    iin=4*(numpt);
    result.write((char*)&iin, sizeof (int));
	PARTLOOP
	if(PP.Flag[n]>=print_flag)
	{
		ffn=float(PP.d50/2);
		result.write((char*)&ffn, sizeof (float));
	}
	

	//  XYZ
	iin=4*(numpt)*3;
	result.write((char*)&iin, sizeof (int));
    PARTLOOP
    if(PP.Flag[n]>=print_flag)
	{
	ffn=float(PP.X[n]);
	result.write((char*)&ffn, sizeof (float));

	ffn=float(PP.Y[n]);
	result.write((char*)&ffn, sizeof (float));

	ffn=float(PP.Z[n]);
	result.write((char*)&ffn, sizeof (float));
	}
	
	//  Connectivity
	count=0;
    iin=4*(numpt);
    result.write((char*)&iin, sizeof (int));
	PARTLOOP
	if(PP.Flag[n]>=print_flag)
	{
	iin=int(count);
	result.write((char*)&iin, sizeof (int));
	++count;
	}

	//  Offset of Connectivity
	count=1;
    iin=4*(numpt);
    result.write((char*)&iin, sizeof (int));
	PARTLOOP
    if(PP.Flag[n]>=print_flag)
	{
	iin=int(count);
	result.write((char*)&iin, sizeof (int));
	++count;
	}

	result<<endl<<"</AppendedData>"<<endl;
    result<<"</VTKFile>"<<endl;

	result.close();
	}

