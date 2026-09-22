/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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
Author: Alexander Hanke
--------------------------------------------------------------------*/

#include"print_gage_location.h"
#include"lexer.h"
#include<sys/stat.h>
#include<sys/types.h>

void print_gage_location::print_wsf_gage_location(lexer *p)
{
    char name[100];
    int iin,offset[100];
    float ffn;
    int count,n;
    int Np = p->P51;
    
    double *Fx = p->P51_x;
	double *Fy = p->P51_y;
    
    double Zloc = p->F60 + 1.1*MAX(p->B91_1,p->B93_1);
    
	mkdir("./REEF3D_Log-Probes",0777);
    

    sprintf(name,"./REEF3D_Log-Probes/REEF3D_wsf_gage_location.vtu");
    
    



	ofstream result;
	result.open(name, ios::binary);

    n=0;

	offset[n]=0;
	++n;
	
	offset[n]=offset[n-1]+4*(Np)+4;
	++n;
	offset[n]=offset[n-1]+4*(Np)+4;
	++n;	
	offset[n]=offset[n-1]+4*(Np)+4;
	++n;	
    offset[n]=offset[n-1]+4*(Np)+4;
	++n;	
	
	// end scalars
    offset[n]=offset[n-1]+4*(Np)*3+4;
    ++n;
    offset[n]=offset[n-1]+4*(Np)*2+4;
    ++n;
	offset[n]=offset[n-1]+4*(Np)+4;
    ++n;
	offset[n]=offset[n-1]+4*(Np)+4;
    ++n;

	//---------------------------------------------
	n=0;
	result<<"<?xml version=\"1.0\"?>"<<endl;
	result<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"<<endl;
	result<<"<UnstructuredGrid>"<<endl;
	result<<"<Piece NumberOfPoints=\""<<Np<<"\" NumberOfCells=\""<<Np<<"\">"<<endl;
	
	
	result<<"<PointData >"<<endl;
    result<<"<DataArray type=\"Float32\" Name=\"radius\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	result<<"<DataArray type=\"Float32\" Name=\"X_Coord\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"Y_Coord\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"Z_Coord\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
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
	

//  radius
    iin=4*(Np);
    result.write((char*)&iin, sizeof (int));
	for(n=0;n<Np;++n)
	{
	ffn=0.1*p->DXM;
	result.write((char*)&ffn, sizeof (float));
	}
	
//  X_coord
    iin=4*(Np);
    result.write((char*)&iin, sizeof (int));
	for(n=0;n<Np;++n)
	{
	ffn=float(p->Xout(Fx[n],Fy[n]));
	result.write((char*)&ffn, sizeof (float));
	}
    
//  Y_Coord
    iin=4*(Np);
    result.write((char*)&iin, sizeof (int));
	for(n=0;n<Np;++n)
	{
	ffn=float(p->Yout(Fx[n],Fy[n]));
	result.write((char*)&ffn, sizeof (float));
	}
    
//  Z_Coord
    iin=4*(Np);
    result.write((char*)&iin, sizeof (int));
	for(n=0;n<Np;++n)
	{
	ffn=float(Zloc);
	result.write((char*)&ffn, sizeof (float));
	}


//  XYZ
	iin=4*(Np)*3;
	result.write((char*)&iin, sizeof (int));
    for(n=0;n<Np;++n)
	{
	ffn=float(p->Xout(Fx[n],Fy[n]));
	result.write((char*)&ffn, sizeof (float));

	ffn=float(p->Yout(Fx[n],Fy[n]));
	result.write((char*)&ffn, sizeof (float));

	ffn=float(Zloc);
	result.write((char*)&ffn, sizeof (float));
	}
	
//  Connectivity
	count=0;
    iin=4*(Np)*2;
    result.write((char*)&iin, sizeof (int));
	for(n=0;n<Np;++n)
	{
	iin=int(0);
	result.write((char*)&iin, sizeof (int));

	iin=int(count);
	result.write((char*)&iin, sizeof (int));
	++count;
	}

//  Offset of Connectivity
	count=0;
    iin=4*(Np);
    result.write((char*)&iin, sizeof (int));
	for(n=0;n<Np;++n)
	{
	iin=(count+1)*2;
	result.write((char*)&iin, sizeof (int));
	++count;
	}


//  Cell types
    iin=4*(Np);
    result.write((char*)&iin, sizeof (int));
	for(n=0;n<Np;++n)
	{
	iin=1;
	result.write((char*)&iin, sizeof (int));
	}

	result<<endl<<"</AppendedData>"<<endl;
    result<<"</VTKFile>"<<endl;

	result.close();
}
