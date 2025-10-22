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
Author: Alexander Hanke
--------------------------------------------------------------------*/

#include "vtr3D.h"
#include "lexer.h"
#include "fdm.h"
#include "ghostcell.h"
#include <sys/stat.h>

void vtr3D::folder(const char* A10)
{
    snprintf(pname,sizeof(pname),"./REEF3D_%s_VTR", A10);
    mkdir(pname,0777);
}

void vtr3D::fileName(char *name, const unsigned int size, const char *A10, const int num, const int rank)
{
    snprintf(name,size,"./REEF3D_%s_VTR/REEF3D-%s-%08i-%06i.vtr",A10,A10,num,rank);
}

void vtr3D::parallelFileName(char *name, const unsigned int size, const char *A10, const int num)
{
    snprintf(name,size,"./REEF3D_%s_VTR/REEF3D-%s-%08i.pvtr",A10,A10,num);
}

void vtr3D::extent(lexer *p, ghostcell *pgc)
{
    int iextent[6];
    iextent[0]=p->origin_i;
    iextent[1]=p->origin_i+p->knox;
    iextent[2]=p->origin_j;
    iextent[3]=p->origin_j+p->knoy;
    iextent[4]=p->origin_k;
    iextent[5]=p->origin_k+p->knoz;

    if(p->mpirank == 0)
        piextent = (int *)malloc(p->mpi_size*6*sizeof(int));
    pgc->gather_int(iextent,6,piextent,6);
}

void vtr3D::offset(lexer *p, int *offset, int &n)
{
    //x
    offset[n]=offset[n-1]+sizeof(float)*(p->knox+1)+sizeof(int); 
    ++n;
    //y
    offset[n]=offset[n-1]+sizeof(float)*(p->knoy+1)+sizeof(int); 
    ++n;
    //z
    offset[n]=offset[n-1]+sizeof(float)*(p->knoz+1)+sizeof(int); 
    ++n;
}

void vtr3D::beginning(lexer *p, std::ostream &result)
{
    xmlVersion(result);
    result<<"<VTKFile type=\"RectilinearGrid\" ";
    vtkVersion(result);
    result<<"<RectilinearGrid WholeExtent=\""<<p->origin_i<<" "<<p->origin_i+p->knox<<" "<<p->origin_j<<" "<<p->origin_j+p->knoy<<" "<<p->origin_k<<" "<<p->origin_k+p->knoz<<" "<<"\">\n";
    if(p->P16==1)
        timeValue(result,p->P34<0.0?p->simtime:p->sedtime);
    result<<"<Piece Extent=\""<<p->origin_i<<" "<<p->origin_i+p->knox<<" "<<p->origin_j<<" "<<p->origin_j+p->knoy<<" "<<p->origin_k<<" "<<p->origin_k+p->knoz<<"\">\n";
}

void vtr3D::beginningParallel(lexer *p, std::ostream &result)
{
    xmlVersion(result);
    result<<"<VTKFile type=\"PRectilinearGrid\" ";
    vtkVersion(result);
    result<<"<PRectilinearGrid WholeExtent=\"0 "<<p->gknox<<" 0 "<<p->gknoy<<" 0 "<<p->gknoz<<"\" GhostLevel=\"0\" Origin=\"0 0 0\" Spacing=\"1 1 1\">\n";
    if(p->P16==1)
        timeValue(result,p->P34<0.0?p->simtime:p->sedtime);
}

void vtr3D::ending(std::ostream &result, const int *offset, int &n)
{
    result<<"<Coordinates>\n";
    result<<"<DataArray type=\"Float32\" Name=\"X\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    n++;
    result<<"<DataArray type=\"Float32\" Name=\"Y\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    n++;
    result<<"<DataArray type=\"Float32\" Name=\"Z\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    n++;
    result<<"</Coordinates>\n";
    result<<"</Piece>\n";
    result<<"</RectilinearGrid>\n";
    appendData(result);
}

void vtr3D::endingParallel(std::ostream &result, const char *A10, const int M10, const int num)
{
    result<<"<PCoordinates>\n";
    result<<"<PDataArray type=\"Float32\" Name=\"X\" format=\"appended\"/>\n";
    result<<"<PDataArray type=\"Float32\" Name=\"Y\" format=\"appended\"/>\n";
    result<<"<PDataArray type=\"Float32\" Name=\"Z\" format=\"appended\"/>\n";
    result<<"</PCoordinates>\n";

    for(int n=0; n<M10; ++n)
    {
        snprintf(pname,sizeof(pname),"REEF3D-%s-%08i-%06i.vtr",A10,num,n+1);
        snprintf(pextent,sizeof(pextent),"%i %i %i %i %i %i",piextent[0+6*n],piextent[1+6*n],piextent[2+6*n],piextent[3+6*n],piextent[4+6*n],piextent[5+6*n]);
        result<<"<Piece Extent=\""<<pextent<<"\" Source=\""<<pname<<"\"/>\n";
    }

    result<<"</PRectilinearGrid>\n";
    result<<"</VTKFile>"<<flush;
}

void vtr3D::structureWrite(lexer *p, fdm*, std::ostream &result)
{
    float ffn;
    int iin;

    // Coordinates
    // x
    iin=sizeof(float)*(p->knox+1);
    result.write((char*)&iin, sizeof (int));
    ITLOOP
    {
        ffn=float(p->XN[IP]);
        result.write((char*)&ffn, sizeof (float));
    }
    // y
    iin=sizeof(float)*(p->knoy+1);
    result.write((char*)&iin, sizeof (int));
    JTLOOP
    {
        ffn=float(p->YN[JP]);
        result.write((char*)&ffn, sizeof (float));
    }
    // z
    iin=sizeof(float)*(p->knoz+1);
    result.write((char*)&iin, sizeof (int));
    KTLOOP
    {
        ffn=float(p->ZN[KP]);
        result.write((char*)&ffn, sizeof (float));
    }

    result<<"\n</AppendedData>\n";
    result<<"</VTKFile>"<<flush;
}