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

#include "vts3D.h"
#include "lexer.h"
#include "fdm.h"
#include "ghostcell.h"
#include <sys/stat.h>

vts3D::vts3D()
{
}

vts3D::~vts3D()
{
}

void vts3D::folder()
{
	mkdir("./REEF3D_CFD_VTS",0777);
}

void vts3D::offset(lexer *p, int *offset, int &n)
{
	// Points
    offset[n]=offset[n-1]+4*(p->pointnum)*3+4;
    ++n;
}

void vts3D::beginning(lexer *p, std::ofstream &result)
{
	xmlVersion(result);
	result<<"<VTKFile type=\"StructuredGrid\" ";
	vtkVersion(result);
	result<<"<StructuredGrid WholeExtent=\""<<p->origin_i<<" "<<p->origin_i+p->knox<<" "<<p->origin_j<<" "<<p->origin_j+p->knoy<<" "<<p->origin_k<<" "<<p->origin_k+p->knoz<<"\">"<<endl;
	result<<"<Piece Extent=\""<<p->origin_i<<" "<<p->origin_i+p->knox<<" "<<p->origin_j<<" "<<p->origin_j+p->knoy<<" "<<p->origin_k<<" "<<p->origin_k+p->knoz<<"\">"<<endl;
}

void vts3D::beginningParallel(lexer *p, std::ofstream &result)
{
	xmlVersion(result);
	result<<"<VTKFile type=\"PStructuredGrid\" ";
	vtkVersion(result);
	result<<"<PStructuredGrid WholeExtent=\"0 "<<p->gknox<<" 0 "<<p->gknoy<<" 0 "<<p->gknoz<<"\" GhostLevel=\"0\" Origin=\"0 0 0\" Spacing=\"1 1 1\">"<<endl;
}

void vts3D::ending(std::ofstream &result, int *offset, int &n)
{
    result<<"<Points>"<<endl;
	result<<"\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset[n]<<"\"/>"<<endl;
	n++;
	result<<"</Points>"<<endl;
	result<<"</Piece>"<<endl;
	result<<"</StructuredGrid>"<<endl;
}

void vts3D::endingParallel(std::ofstream &result, int &M10, int &num)
{
	result<<"<PPoints>"<<endl;
	result<<"\t<PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>"<<endl;
    result<<"</PPoints>"<<endl;

	for(int n=0; n<M10; ++n)
	{
		sprintf(pname,"REEF3D-CFD-%08i-%06i.vts",num,n+1);
		sprintf(pextent,"%i %i %i %i %i %i",piextent[0+6*n],piextent[1+6*n],piextent[2+6*n],piextent[3+6*n],piextent[4+6*n],piextent[5+6*n]);
		result<<"<Piece Extent=\""<<pextent<<"\" Source=\""<<pname<<"\"/>"<<endl;
	}

	result<<"</PStructuredGrid>"<<endl;
	result<<"</VTKFile>"<<flush;
}

void vts3D::extent(lexer *p, ghostcell *pgc)
{
	int iextent[6];
	iextent[0]=p->origin_i;
	iextent[1]=p->origin_i+p->knox;
	iextent[2]=p->origin_j;
	iextent[3]=p->origin_j+p->knoy;
	iextent[4]=p->origin_k;
	iextent[5]=p->origin_k+p->knoz;
	
    if ( p->mpirank == 0)
    piextent = (int *)malloc(p->mpi_size*6*sizeof(int)); 
    pgc->gather_int(iextent,6,piextent,6);
}

void vts3D::structureWrite(lexer *p, fdm *a, std::ofstream &result)
{
	float ffn;
	int iin;

	// Coordinates
	// points
	iin=3*4*(p->pointnum);
	result.write((char*)&iin, sizeof (int));
    KTLOOP
	JTLOOP
	ITLOOP
	TPCHECK
	{
	ffn=float(p->XN[IP]);
	result.write((char*)&ffn, sizeof (float));

	ffn=float(p->YN[JP]);
	result.write((char*)&ffn, sizeof (float));

	ffn=float(p->ZN[KP]);
	result.write((char*)&ffn, sizeof (float));
	}

	result<<endl<<"</AppendedData>"<<endl;
    result<<"</VTKFile>"<<flush;
}