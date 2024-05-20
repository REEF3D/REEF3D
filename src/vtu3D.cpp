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

#include "vtu3D.h"
#include "lexer.h"
#include "fdm.h"
#include <sys/stat.h>

vtu3D::vtu3D()
{
}

vtu3D::~vtu3D()
{
}

void vtu3D::folder()
{
	mkdir("./REEF3D_CFD_VTU",0777);
}

void vtu3D::offset(lexer *p, int *offset, int &n)
{
	// Points
    offset[n]=offset[n-1]+4*(p->pointnum)*3+4;
    ++n;

	// Cells
    offset[n]=offset[n-1] + 4*(p->tpcellnum*8) + 4;
    ++n;
    offset[n]=offset[n-1] + 4*(p->tpcellnum) + 4;
    ++n;
	offset[n]=offset[n-1] + 4*(p->tpcellnum) + 4;
    ++n;
}

void vtu3D::beginning(lexer *p, std::ofstream &result)
{
	xmlVersion(result);
	result<<"<VTKFile type=\"UnstructuredGrid\" ";
	vtkVersion(result);
	result<<"<UnstructuredGrid>"<<endl;
	result<<"<Piece NumberOfPoints=\""<<p->pointnum<<"\" NumberOfCells=\""<<p->tpcellnum<<"\">"<<endl;
}

void vtu3D::beginningParallel(lexer *p, std::ofstream &result)
{
	xmlVersion(result);
	result<<"<VTKFile type=\"PUnstructuredGrid\" ";
	vtkVersion(result);
	result<<"<PUnstructuredGrid GhostLevel=\"0\">"<<endl;
}

void vtu3D::ending(std::ofstream &result, int *offset, int &n)
{
    result<<"<Points>"<<endl;
    result<<"\t<DataArray type=\"Float32\"  NumberOfComponents=\"3\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"</Points>"<<endl;

    result<<"<Cells>"<<endl;
    result<<"\t<DataArray type=\"Int32\"  Name=\"connectivity\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	result<<"\t<DataArray type=\"Int32\"  Name=\"offsets\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
	++n;
    result<<"\t<DataArray type=\"Int32\"  Name=\"types\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	result<<"</Cells>"<<endl;

    result<<"</Piece>"<<endl;
    result<<"</UnstructuredGrid>"<<endl;
}

void vtu3D::endingParallel(std::ofstream &result, int &M10, int &num)
{
	result<<"<PPoints>"<<endl;
	result<<"\t<PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>"<<endl;
	result<<"</PPoints>"<<endl;

	result<<"<Cells>"<<endl;
	result<<"\t<DataArray type=\"Int32\"  Name=\"connectivity\"/>"<<endl;
	result<<"\t<DataArray type=\"Int32\"  Name=\"offsets\" />"<<endl;
	result<<"\t<DataArray type=\"Int32\"  Name=\"types\" />"<<endl;
	result<<"</Cells>"<<endl;

	char pname[30];
	for(int n=0; n<M10; ++n)
	{
	sprintf(pname,"REEF3D-CFD-%08i-%06i.vtu",num,n+1);
	result<<"<Piece Source=\""<<pname<<"\"/>"<<endl;
	}

	result<<"</PUnstructuredGrid>"<<endl;
	result<<"</VTKFile>"<<flush;
}

void vtu3D::structureWrite(lexer *p, fdm *a, std::ofstream &result)
{
	float ffn;
	int iin;
	double phase=0.0;
	double zcoor;

	//  XYZ
	double theta_y = p->B192_1*(PI/180.0);
	double omega_y = 2.0*PI*p->B192_2;

	if(p->B192==1 && p->simtime>=p->B194_s && p->simtime<=p->B194_e)
	phase = omega_y*p->simtime;

	iin=4*(p->pointnum)*3;
	result.write((char*)&iin, sizeof (int));
	TPLOOP
	{

		zcoor=p->ZN[KP1];

		ffn=float( (p->XN[IP1]-p->B192_3)*cos(theta_y*sin(phase)) - (zcoor-p->B192_4)*sin(theta_y*sin(phase)) + p->B192_3 
					+ p->B181_1*sin((2.0*PI*p->B181_2)*p->simtime + p->B181_3));
		result.write((char*)&ffn, sizeof (float));

		ffn=float(p->YN[JP1]) + p->B182_1*std::sin((2.0*PI*p->B182_2)*p->simtime + p->B182_3);
		result.write((char*)&ffn, sizeof (float));
		

		ffn=float((p->XN[IP1]-p->B192_3)*sin(theta_y*sin(phase)) + (zcoor-p->B192_4)*cos(theta_y*sin(phase)) + p->B192_4
					+ p->B183_1*sin((2.0*PI*p->B183_2)*p->simtime + p->B183_3));
		result.write((char*)&ffn, sizeof (float));
	}

	//  Connectivity
    iin=4*(p->tpcellnum)*8;
    result.write((char*)&iin, sizeof (int));
    BASEREVLOOP
    if(p->flag5[IJK]!=-20 && p->flag5[IJK]!=-30)
	{
		iin=int(a->nodeval(i-1,j-1,k-1)-1);
		result.write((char*)&iin, sizeof (int));

		iin=int(a->nodeval(i,j-1,k-1))-1;
		result.write((char*)&iin, sizeof (int));

		iin= int(a->nodeval(i,j,k-1))-1;
		result.write((char*)&iin, sizeof (int));

		iin=int(a->nodeval(i-1,j,k-1))-1;
		result.write((char*)&iin, sizeof (int));

		iin=int(a->nodeval(i-1,j-1,k))-1;
		result.write((char*)&iin, sizeof (int));

		iin=int(a->nodeval(i,j-1,k))-1;
		result.write((char*)&iin, sizeof (int));

		iin=int(a->nodeval(i,j,k))-1;
		result.write((char*)&iin, sizeof (int));

		iin=int(a->nodeval(i-1,j,k))-1;
		result.write((char*)&iin, sizeof (int));
	}

	// //  Offset of Connectivity
    iin=4*(p->tpcellnum);
    result.write((char*)&iin, sizeof (int));
	for(n=0;n<p->tpcellnum;++n)
	{
		iin=(n+1)*8;
		result.write((char*)&iin, sizeof (int));
	}

	//  Cell types
    iin=4*(p->tpcellnum);
    result.write((char*)&iin, sizeof (int));
	for(n=0;n<p->tpcellnum;++n)
	{
		iin=12;
		result.write((char*)&iin, sizeof (int));
	}

	result<<endl;
	result<<"</AppendedData>"<<endl;
    result<<"</VTKFile>"<<flush;
}