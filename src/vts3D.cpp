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
#include "fdm_fnpf.h"
#include "fdm_nhf.h"
#include "ghostcell.h"
#include <sys/stat.h>

vts3D::vts3D()
{
}

vts3D::~vts3D()
{
}

void vts3D::folder(const char *A10)
{
	sprintf(pname,"./REEF3D_%s_VTS", A10);
	mkdir(pname,0777);
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
	result<<"<StructuredGrid WholeExtent=\""<<p->origin_i<<" "<<p->origin_i+p->knox<<" "<<p->origin_j<<" "<<p->origin_j+p->knoy<<" "<<p->origin_k<<" "<<p->origin_k+p->knoz<<"\">\n";
	if(p->P16==1)
    {
        result<<"<FieldData>\n";
        result<<"<DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\"> "<<p->simtime<<"\n";
        result<<"</DataArray>\n";
        result<<"</FieldData>\n";
    }
    result<<"<Piece Extent=\""<<p->origin_i<<" "<<p->origin_i+p->knox<<" "<<p->origin_j<<" "<<p->origin_j+p->knoy<<" "<<p->origin_k<<" "<<p->origin_k+p->knoz<<"\">\n";
}

void vts3D::beginningParallel(lexer *p, std::ofstream &result)
{
	xmlVersion(result);
	result<<"<VTKFile type=\"PStructuredGrid\" ";
	vtkVersion(result);
	result<<"<PStructuredGrid WholeExtent=\"0 "<<p->gknox<<" 0 "<<p->gknoy<<" 0 "<<p->gknoz<<"\" GhostLevel=\"0\" Origin=\"0 0 0\" Spacing=\"1 1 1\">\n";
}

void vts3D::ending(std::ofstream &result, const int *offset, int &n)
{
    result<<"<Points>\n";
	result<<"\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
	n++;
	result<<"</Points>\n";
	result<<"</Piece>\n";
	result<<"</StructuredGrid>\n";
}

void vts3D::endingParallel(std::ofstream &result, const char *A10, const int M10, const int num)
{
	result<<"<PPoints>\n";
	result<<"\t<PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>\n";
    result<<"</PPoints>\n";

	for(int n=0; n<M10; ++n)
	{
		sprintf(pname,"REEF3D-%s-%08i-%06i.vts",A10,num,n+1);
		sprintf(pextent,"%i %i %i %i %i %i",piextent[0+6*n],piextent[1+6*n],piextent[2+6*n],piextent[3+6*n],piextent[4+6*n],piextent[5+6*n]);
		result<<"<Piece Extent=\""<<pextent<<"\" Source=\""<<pname<<"\"/>\n";
	}

	result<<"</PStructuredGrid>\n";
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

	result<<"\n"<<"</AppendedData>\n";
    result<<"</VTKFile>"<<flush;
}

void vts3D::structureWrite(lexer *p, fdm_fnpf *c, std::ofstream &result)
{
	float ffn;
	int iin;
	double phase;
	double zcoor;

	//  XYZ
	iin=4*(p->pointnum)*3;
	result.write((char*)&iin, sizeof (int));

    double theta_y = p->B192_1*(PI/180.0);
	double omega_y = 2.0*PI*p->B192_2;
    double waterlevel;

    if(p->B192==1 && p->simtime>=p->B194_s && p->simtime<=p->B194_e)
    phase = omega_y*p->simtime;

    TPLOOP
	{
        waterlevel = p->sl_ipol4eta(p->wet,c->eta,c->bed)+p->wd - p->sl_ipol4(c->bed);

        zcoor = p->ZN[KP1]*waterlevel + p->sl_ipol4(c->bed);

        if(p->wet[IJ]==0)
        zcoor=c->bed(i,j);

        if(i+p->origin_i==-1 && j+p->origin_j==-1 && p->wet[(0-p->imin)*p->jmax + (0-p->jmin)]==1)
        zcoor = p->ZN[KP1]*c->WL(i,j) + c->bed(i,j);

        

        ffn=float( (p->XN[IP1]-p->B192_3)*cos(theta_y*sin(phase)) - (zcoor-p->B192_4)*sin(theta_y*sin(phase)) + p->B192_3);
        result.write((char*)&ffn, sizeof (float));

        ffn=float(p->YN[JP1]);
        result.write((char*)&ffn, sizeof (float));

        ffn=float((p->XN[IP1]-p->B192_3)*sin(theta_y*sin(phase)) + (zcoor-p->B192_4)*cos(theta_y*sin(phase)) + p->B192_4);
        result.write((char*)&ffn, sizeof (float));
	}

	result<<"\n"<<"</AppendedData>\n";
    result<<"</VTKFile>"<<flush;
}

void vts3D::structureWrite(lexer *p, fdm_nhf *d, std::ofstream &result)
{
	float ffn;
	int iin;
	double phase;
	double zcoor;

	//  XYZ
	iin=4*(p->pointnum)*3;
	result.write((char*)&iin, sizeof (int));
    TPLOOP
	{
        zcoor = p->ZN[KP1]*p->sl_ipol4(d->WL) + p->sl_ipol4(d->bed);

        if(p->wet[IJ]==0)
        zcoor=p->sl_ipol4(d->bed);
        
        if(i+p->origin_i==-1 && j+p->origin_j==-1 && p->wet[(0-p->imin)*p->jmax + (0-p->jmin)]==1)
        zcoor = p->ZN[KP1]*d->WL(i,j) + d->bed(i,j);

        // -- 
        ffn=float(p->XN[IP1]);
        result.write((char*)&ffn, sizeof (float));

        ffn=float(p->YN[JP1]);
        result.write((char*)&ffn, sizeof (float));

        ffn=float(zcoor);
        result.write((char*)&ffn, sizeof (float));
	}

	result<<"\n"<<"</AppendedData>\n";
    result<<"</VTKFile>"<<flush;
}