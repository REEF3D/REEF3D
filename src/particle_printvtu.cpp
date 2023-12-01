/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"particle_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void particle_f::print_particles(lexer* p, fdm* a, ghostcell* pgc)
{
    if(((p->count%p->Q181==0 && p->Q182<0.0 && p->Q180==1 )|| (p->count==0 &&  p->Q182<0.0 && p->Q180==1)) && p->Q181>0)
	{
    print_vtu(p,a,pgc,pos,posflag,posactive,1);
	++printcount;
	}
    
    if((p->simtime>p->fsfprinttime && p->Q182>0.0 && p->Q180==1) || (p->count==0 &&  p->Q182>0.0))
    {
    print_vtu(p,a,pgc,pos,posflag,posactive,1);
    p->partprinttime+=p->Q182;
    }
    
}

void particle_f::print_vtu(lexer* p, fdm* a, ghostcell* pgc,double** f,int *flag,int active, int sign)
{
	int numpt=0;
	int count;
	
	for(n=0;n<active;++n)
    if(flag[n]>0)
	++numpt;
	
    cout<<"PACTIVE-"<<p->mpirank<<": "<<numpt<<" "<<active<<endl;
	
	if(p->mpirank==0)
	pvtu_pos(a,p,pgc);


    header_pos(a,p,pgc);

	ofstream result;
	result.open(name, ios::binary);

    n=0;

	offset[n]=0;
	++n;
	
	offset[n]=offset[n-1]+4*(numpt)+4;
	++n;
	offset[n]=offset[n-1]+4*(numpt)+4;
	++n;	
	offset[n]=offset[n-1]+4*(numpt)+4;
	++n;	
	
	// end scalars
    offset[n]=offset[n-1]+4*(numpt)*3+4;
    ++n;
    offset[n]=offset[n-1]+4*(numpt)*2+4;
    ++n;
	offset[n]=offset[n-1]+4*(numpt)+4;
    ++n;
	offset[n]=offset[n-1]+4*(numpt)+4;
    ++n;

	//---------------------------------------------
	n=0;
	result<<"<?xml version=\"1.0\"?>"<<endl;
	result<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"<<endl;
	result<<"<UnstructuredGrid>"<<endl;
	result<<"<Piece NumberOfPoints=\""<<numpt<<"\" NumberOfCells=\""<<numpt<<"\">"<<endl;
	
	
	result<<"<PointData >"<<endl;
	result<<"<DataArray type=\"Float32\" Name=\"phi\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"radius\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	result<<"<DataArray type=\"Float32\" Name=\"correction\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	result<<"</PointData>"<<endl;
	
	

    result<<"<Points>"<<endl;
    result<<"<DataArray type=\"Float32\"  NumberOfComponents=\"3\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"</Points>"<<endl;
	
	

    result<<"<Cells>"<<endl;
	result<<"<DataArray type=\"Int32\"  Name=\"connectivity\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	result<<"<DataArray type=\"Int32\"  Name=\"offsets\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
	++n;
    result<<"<DataArray type=\"Int32\"  Name=\"types\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	result<<"</Cells>"<<endl;

    result<<"</Piece>"<<endl;
    result<<"</UnstructuredGrid>"<<endl;

//----------------------------------------------------------------------------
    result<<"<AppendedData encoding=\"raw\">"<<endl<<"_";
	

//  lsv
    iin=4*(numpt);
    result.write((char*)&iin, sizeof (int));
	for(n=0;n<active;++n)
    if(flag[n]>0)
	{
	ffn=float(f[n][3]);
	result.write((char*)&ffn, sizeof (float));
	}
	
//  radius
    iin=4*(numpt);
    result.write((char*)&iin, sizeof (int));
	for(n=0;n<active;++n)
    if(flag[n]>0)
	{
	ffn=float(f[n][4]);
	result.write((char*)&ffn, sizeof (float));
	}

//  correction
    iin=4*(numpt);
    result.write((char*)&iin, sizeof (int));
	for(n=0;n<active;++n)
    if(flag[n]>0)
	{
		if(sign==1)
		{
		if(f[n][3]<=-f[n][4])
		ffn=float(1.0);
		
		if(f[n][3]>-f[n][4])
		ffn=float(0.0);
		}
		
		if(sign==2)
		{
		if(f[n][3]>=f[n][4])
		ffn=float(1.0);
		
		if(f[n][3]<f[n][4])
		ffn=float(0.0);
		}
		
	result.write((char*)&ffn, sizeof (float));
	}

//  XYZ
	iin=4*(numpt)*3;
	result.write((char*)&iin, sizeof (int));
    for(n=0;n<active;++n)
    if(flag[n]>0)
	{
	ffn=float(f[n][0]);
	result.write((char*)&ffn, sizeof (float));

	ffn=float(f[n][1]);
	result.write((char*)&ffn, sizeof (float));

	ffn=float(f[n][2]);
	result.write((char*)&ffn, sizeof (float));
	}
	
//  Connectivity
	count=0;
    iin=4*(numpt)*2;
    result.write((char*)&iin, sizeof (int));
	for(n=0;n<active;++n)
	if(flag[n]>0)
	{
	iin=int(0);
	result.write((char*)&iin, sizeof (int));

	iin=int(count);
	result.write((char*)&iin, sizeof (int));
	++count;
	}

//  Offset of Connectivity
	count=0;
    iin=4*(numpt);
    result.write((char*)&iin, sizeof (int));
	for(n=0;n<active;++n)
    if(flag[n]>0)
	{
	iin=(count+1)*2;
	result.write((char*)&iin, sizeof (int));
	++count;
	}


//  Cell types
    iin=4*(numpt);
    result.write((char*)&iin, sizeof (int));
	for(n=0;n<active;++n)
    if(flag[n]>0)
	{
	iin=1;
	result.write((char*)&iin, sizeof (int));
	}

	result<<endl<<"</AppendedData>"<<endl;
    result<<"</VTKFile>"<<endl;

	result.close();
}

