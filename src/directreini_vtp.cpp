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

#include"directreini.h"
#include"lexer.h"

void directreini::vtp(lexer* p)
{
    int num=0;
    if(p->P15==1)
        num = p->printcount;
    else if(p->P15==2)
        num = p->count;

    if(p->mpirank==0)
        pvtp(p,num);

    sprintf(name,"./REEF3D_FSF/REEF3D-FSF-%08i-%06i.vtp",num,p->mpirank+1);

    ofstream result;
    result.open(name, ios::binary);
    //---------------------------------------------

    polygon_num=facount;

    polygon_sum=0;
    for(n=0;n<polygon_num;++n)
        polygon_sum+=numfac[n];

    vertice_num = ccptcount;

    cout<<p->mpirank<<" Vertice_num: "<<vertice_num<<" Polygon_Num: "<<polygon_num<<" Polygon_Sum: "<<polygon_sum<<endl;

    //---------------------------------------------
    n=0;
    offset[n]=0;
    ++n;

    offset[n]=offset[n-1] + sizeof(float)*vertice_num*3 + sizeof(int);
    ++n;
    offset[n]=offset[n-1] + sizeof(int)*polygon_sum + sizeof(int);
    ++n;
    offset[n]=offset[n-1] + sizeof(int)*polygon_num + sizeof(int);
    ++n;
    //---------------------------------------------

    result<<"<?xml version=\"1.0\"?>\n";
    result<<"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    result<<"<PolyData>\n";
    result<<"<Piece NumberOfPoints=\""<<vertice_num<<"\" NumberOfPolys=\""<<polygon_num<<"\">\n";

    n=0;
    result<<"<Points>\n";
    result<<"<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"</Points>\n";

    result<<"<Polys>\n";
    result<<"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"</Polys>\n";

    result<<"</Piece>\n";
    result<<"</PolyData>\n";
    result<<"<AppendedData encoding=\"raw\">\n_";

    //----------------------------------------------------------------------------

    float ffn;
    int iin;
    //  XYZ
    iin=4*vertice_num*3;
    result.write((char*)&iin, sizeof (int));
    for(n=0;n<vertice_num;++n)
    {
        ffn=ccpt[n][0];
        result.write((char*)&ffn, sizeof (float));

        ffn=ccpt[n][1];
        result.write((char*)&ffn, sizeof (float));

        ffn=ccpt[n][2];
        result.write((char*)&ffn, sizeof (float));
    }

    //  Connectivity POLYGON
    iin=4*polygon_sum;
    result.write((char*)&iin, sizeof (int));
    for(n=0;n<polygon_num;++n)
    for(q=0;q<numfac[n];++q)
    {
        iin=facet[n][q];
        result.write((char*)&iin, sizeof (int));
    }

    //  Offset of Connectivity
    iin=4*polygon_num;
    result.write((char*)&iin, sizeof (int));
    iin=0;
    for(n=0;n<polygon_num;++n)
    {
        iin+= numfac[n];
        result.write((char*)&iin, sizeof (int));
    }

    result<<"\n</AppendedData>\n";
    result<<"</VTKFile>\n";

    result.close();
}
