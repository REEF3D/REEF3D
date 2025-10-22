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

#include"force.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void force::print_vtp(lexer* p, fdm* a, ghostcell *pgc)
{
    polygon_num=facount;

    polygon_sum=0;
    for(n=0;n<polygon_num;++n)
        polygon_sum+=numpt[n];

    vertice_num = ccptcount;

    //---------------------------------------------
    n=0;
    offset[n]=0;
    ++n;
    offset[n]=offset[n-1] + sizeof(float)*vertice_num*3 + sizeof(int);
    ++n;
    //Data
    offset[n]=offset[n-1] + sizeof(float)*vertice_num*3 + sizeof(int);
    ++n;
    offset[n]=offset[n-1] + sizeof(float)*vertice_num + sizeof(int);
    ++n;
    //End Data
    offset[n]=offset[n-1] + sizeof(int)*polygon_sum + sizeof(int);
    ++n;
    offset[n]=offset[n-1] + sizeof(int)*polygon_num + sizeof(int);
    ++n;
    //---------------------------------------------

    int num=0;
    if(p->P15==1)
        num = forceprintcount;
    else if(p->P15==2)
        num = p->count;

    if(p->mpirank==0)
        pvtp(p,num);

    sprintf(name,"./REEF3D_SOLID/REEF3D-SOLID-%i-%08i-%06i.vtp",ID,num,p->mpirank+1);

    ofstream result;
    result.open(name, ios::binary);

    result<<"<?xml version=\"1.0\"?>\n";
    result<<"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    result<<"<PolyData>\n";
    result<<"<Piece NumberOfPoints=\""<<vertice_num<<"\" NumberOfPolys=\""<<polygon_num<<"\">\n";

    n=0;
    result<<"<Points>\n";
    result<<"<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"</Points>\n";

    result<<"<PointData>\n";
    result<<"<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"pressure\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"</PointData>\n";

    result<<"<Polys>\n";
    result<<"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"</Polys>\n";

    result<<"</Piece>\n";
    result<<"</PolyData>\n";
    result<<"<AppendedData encoding=\"raw\">\n"<<"_";

    //----------------------------------------------------------------------------

    //  XYZ
    iin=sizeof(float)*vertice_num*3;
    result.write((char*)&iin, sizeof(int));
    for(n=0;n<vertice_num;++n)
    {
        ffn=ccpt[n][0];
        result.write((char*)&ffn, sizeof(float));

        ffn=ccpt[n][1];
        result.write((char*)&ffn, sizeof(float));

        ffn=ccpt[n][2];
        result.write((char*)&ffn, sizeof(float));
    }

    //  Velocity
    iin=sizeof(float)*vertice_num*3;
    result.write((char*)&iin, sizeof(int));
    for(n=0;n<vertice_num;++n)
    {
        ffn=float(p->ccipol1(a->u,ccpt[n][0],ccpt[n][1],ccpt[n][2]));
        result.write((char*)&ffn, sizeof(float));

        ffn=float(p->ccipol2(a->v,ccpt[n][0],ccpt[n][1],ccpt[n][2]));
        result.write((char*)&ffn, sizeof(float));

        ffn=float(p->ccipol3(a->w,ccpt[n][0],ccpt[n][1],ccpt[n][2]));
        result.write((char*)&ffn, sizeof(float));
    }

    //  Pressure
    iin=sizeof(float)*vertice_num;
    result.write((char*)&iin, sizeof(int));
    for(n=0;n<vertice_num;++n)
    {
        ffn=float(p->ccipol4(a->press,ccpt[n][0],ccpt[n][1],ccpt[n][2]) - p->pressgage);
        result.write((char*)&ffn, sizeof(float));
    }

    //  Connectivity POLYGON
    iin=sizeof(int)*polygon_sum;
    result.write((char*)&iin, sizeof(int));
    for(n=0;n<polygon_num;++n)
    {
        if(numpt[n]==3)
        {
            iin=facet[n][0];
            result.write((char*)&iin, sizeof(int));

            iin=facet[n][1];
            result.write((char*)&iin, sizeof(int));

            iin=facet[n][2];
            result.write((char*)&iin, sizeof(int));
        }
        else if(numpt[n]==4)
        {
            iin=facet[n][0];
            result.write((char*)&iin, sizeof(int));

            iin=facet[n][1];
            result.write((char*)&iin, sizeof(int));

            iin=facet[n][3];
            result.write((char*)&iin, sizeof(int));

            iin=facet[n][2];
            result.write((char*)&iin, sizeof(int));
        }
    }

    //  Offset of Connectivity
    iin=sizeof(int)*polygon_num;
    result.write((char*)&iin, sizeof(int));
    iin=0;
    for(n=0;n<polygon_num;++n)
    {
        iin+= numpt[n];
        result.write((char*)&iin, sizeof(int));
    }

    result<<"\n</AppendedData>\n";
    result<<"</VTKFile>\n";

    result.close();
}
