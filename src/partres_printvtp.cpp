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
for more details->

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#include"partres.h"
#include"lexer.h"
#include"sediment_fdm.h"
#include"ghostcell.h"
#include<fstream>
#include<sstream>
#include<sys/types.h>
#include<cstdio>
#include<cstring>
#include<vector>

void partres::print_particles(lexer* p, sediment_fdm *s)
{
    if((p->count%p->Q181==0 || p->count==0) && (p->Q180==1 && p->Q181>0 && p->Q182<0.0))
    {
        print_vtp(p,s);
        ++printcount;
    }

    if((p->simtime>printtime || p->count==0) && (p->Q180==1 && p->Q181<0 && p->Q182>0.0))
    {
        print_vtp(p,s);
        printtime+=p->Q182;
        ++printcount;
    }
}

void partres::print_vtp(lexer* p, sediment_fdm *s)
{
    int numpt=0;

    for(n=0;n<P.index;++n)
    if(P.Flag[n]>0)
        numpt++;

    int num = printcount;

    if(p->mpirank==0)
        pvtp(p,num);

    int n=0;
    int offset[100];

    offset[n]=0;
    ++n;

    offset[n]=offset[n-1]+numpt*sizeof(float)+sizeof(int); //flag
    ++n;
    if(p->P23==1)
    {
        offset[n]=offset[n-1]+numpt*sizeof(float)+sizeof(int); //Test
        ++n;
    }
    offset[n]=offset[n-1]+3*numpt*sizeof(float)+sizeof(int); //velocity
    ++n;
    offset[n]=offset[n-1]+numpt*sizeof(float)+sizeof(int); //radius
    ++n;
    offset[n]=offset[n-1]+3*numpt*sizeof(float)+sizeof(int); //fluid velocity
    ++n;
    offset[n]=offset[n-1]+numpt*sizeof(float)+sizeof(int); //bedChange
    ++n;

    offset[n]=offset[n-1]+3*numpt*sizeof(float)+sizeof(int); //xyz
    ++n;
    offset[n]=offset[n-1]+numpt*sizeof(int)+sizeof(int); //connectivitey
    ++n;
    offset[n]=offset[n-1]+numpt*sizeof(int)+sizeof(int); //offset connectivity
    ++n;

    //---------------------------------------------

    std::stringstream result;
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

    result<<"<PointData>\n";
    result<<"<DataArray type=\"Float32\" Name=\"Flag\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    if(p->P23==1)
    {
        result<<"<DataArray type=\"Float32\" Name=\"Test\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
        ++n;
    }
    result<<"<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"radius\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"fluid velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"bedChange\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"</PointData>\n";

    vtp3D::points(result,offset,n);

    vtp3D::verts(result,offset,n);

    vtp3D::ending(result);

    //----------------------------------------------------------------------------

    size_t m=result.str().length();

    const size_t total_size = m + offset[n] + 27;

    std::vector<char> buffer(total_size);
    std::memcpy(&buffer[0],result.str().data(),m);
    int iin;
    float ffn;

    // flag
    iin=numpt*sizeof(float);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
    for(n=0;n<P.index;++n)
        if(P.Flag[n]>=0)
        {
            ffn=float(P.Flag[n]);
            std::memcpy(&buffer[m],&ffn,sizeof(float));
            m+=sizeof(float);
        }

    // Test
    if(p->P23==1)
    {
        iin=numpt*sizeof(float);
        std::memcpy(&buffer[m],&iin,sizeof(int));
        m+=sizeof(int);
        for(n=0;n<P.index;++n)
            if(P.Flag[n]>=0)
            {
                ffn=float(P.Test[n]);
                std::memcpy(&buffer[m],&ffn,sizeof(float));
                m+=sizeof(float);
            }
    }

    // velocities
    iin=3*numpt*sizeof(float);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
    for(n=0;n<P.index;++n)
        if(P.Flag[n]>=0)
        {
            ffn=float(P.U[n]);
            std::memcpy(&buffer[m],&ffn,sizeof(float));
            m+=sizeof(float);

            ffn=float(P.V[n]);
            std::memcpy(&buffer[m],&ffn,sizeof(float));
            m+=sizeof(float);

            ffn=float(P.W[n]);
            std::memcpy(&buffer[m],&ffn,sizeof(float));
            m+=sizeof(float);
        }

    // radius
    iin=numpt*sizeof(float);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
    for(n=0;n<P.index;++n)
        if(P.Flag[n]>=0)
        {
            ffn=float(P.d50/2);
            std::memcpy(&buffer[m],&ffn,sizeof(float));
            m+=sizeof(float);
        }

    // fluid velocities
    iin=3*numpt*sizeof(float);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
    for(n=0;n<P.index;++n)
        if(P.Flag[n]>=0)
        {
            ffn=float(P.Uf[n]);
            std::memcpy(&buffer[m],&ffn,sizeof(float));
            m+=sizeof(float);

            ffn=float(P.Vf[n]);
            std::memcpy(&buffer[m],&ffn,sizeof(float));
            m+=sizeof(float);

            ffn=float(P.Wf[n]);
            std::memcpy(&buffer[m],&ffn,sizeof(float));
            m+=sizeof(float);
        }

    // bedChange
    iin=numpt*sizeof(float);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
    for(n=0;n<P.index;++n)
        if(P.Flag[n]>=0)
        {
            //ffn=float(p->ccslipol4(s->bedch,P.X[n],P.Y[n]));
            ffn=float(p->ccslipol4(s->bedzh,P.X[n],P.Y[n])-p->ccslipol4(s->bedzh0,P.X[n],P.Y[n]));
            std::memcpy(&buffer[m],&ffn,sizeof(float));
            m+=sizeof(float);
        }

    //  XYZ
    iin=3*numpt*sizeof(float);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
    for(n=0;n<P.index;++n)
        if(P.Flag[n]>=0)
        {
            ffn=float(P.X[n]);
            std::memcpy(&buffer[m],&ffn,sizeof(float));
            m+=sizeof(float);

            ffn=float(P.Y[n]);
            std::memcpy(&buffer[m],&ffn,sizeof(float));
            m+=sizeof(float);

            ffn=float(P.Z[n]);
            std::memcpy(&buffer[m],&ffn,sizeof(float));
            m+=sizeof(float);
        }

    //  Connectivity
    int count=0;
    iin=numpt*sizeof(int);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
    for(n=0;n<P.index;++n)
        if(P.Flag[n]>=0)
        {
            iin=int(count);
            std::memcpy(&buffer[m],&iin,sizeof(int));
            m+=sizeof(int);
            ++count;
        }

    //  Offset of Connectivity
    count=1;
    iin=numpt*sizeof(int);
    std::memcpy(&buffer[m],&iin,sizeof(int));
    m+=sizeof(int);
    for(n=0;n<P.index;++n)
        if(P.Flag[n]>=0)
        {
            iin=int(count);
            std::memcpy(&buffer[m],&iin,sizeof(int));
            m+=sizeof(int);
            ++count;
        }

    std::stringstream footer;
    vtp3D::footer(footer);
    std::memcpy(&buffer[m],footer.str().data(),footer.str().size());

    char filename[200];
    sprintf(filename,"./REEF3D_CFD_SedPart/REEF3D-SedPart-%08i-%06i.vtp",printcount,p->mpirank+1);

    const size_t SMALL_FILE_THRESHOLD = 10 * 1024 * 1024;   // 10MB

    if(total_size < SMALL_FILE_THRESHOLD)
    {
        // Small files: Simple fwrite with modest buffer
        FILE* file = fopen(filename, "wb");
        if(file)
        {
            setvbuf(file, nullptr, _IOFBF, 32768); // 32KB buffer
            fwrite(buffer.data(), buffer.size(), 1, file);
            fclose(file);
        }
    }
    else
    {
        // Medium files: Optimized C I/O with larger buffer
        FILE* file = fopen(filename, "wb");
        if(file)
        {
            setvbuf(file, nullptr, _IOFBF, 131072); // 128KB buffer
            fwrite(buffer.data(), buffer.size(), 1, file);
            fclose(file);
        }
    }
}
